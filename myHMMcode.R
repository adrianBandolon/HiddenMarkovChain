forward.prob = function (x, A, B, w) {
  # Output the forward probability matrix alp
  # alp: T by mz, (t, i) entry = P(x_{1:t}, Z_t = i)
  
  T = length(x)
  mz = nrow(A)
  alp = matrix(0, T, mz)
  
  # fill in the first row of alp
  alp[1,] = w * B[, x[1]]
  
  # Recursively compute the remaining rows of alp
  for (t in 2:T) {
    tmp = alp[t - 1,] %*% A
    alp[t,] = tmp * B[, x[t]]
  }
  
  return(alp)
}

backward.prob = function (x, A, B, w) {
  # Output the backward probability matrix beta
  # beta: T by mz, (t, i) entry = P(x_{(t+1):n} | Z_t = i)
  # for t=1, ..., n-1
  
  T = length(x)
  mz = nrow(A)
  beta = matrix(1, T, mz)
  
  # The last row of beta is all 1.
  # Recursively compute the previous rows of beta
  for (t in (T - 1):1) {
    tmp = as.matrix(beta[t + 1,] * B[, x[t + 1]])  # make tmp a column vector
    beta[t,] = t(A %*% tmp)
  }
  
  return(beta)
}

BW.onestep = function (x, A, B, w) {
  # Input:
  # x: T-by-1 observation sequence
  # A: current estimate for mz-by-mz transition matrix
  # B: current estimate for mz-by-mx emission matrix
  # w: current estimate for mz-by-1 initial distribution over Z_1
  # Output the updated parameters
  # para = list(A = A1, B = B1)
  
  # We DO NOT update the initial distribution w
  
  T = length(x)
  mz = nrow(A)
  alp = forward.prob(x, A, B, w)
  beta = backward.prob(x, A, B, w)
  myGamma = array(0, dim = c(mz, mz, T - 1))
  
  for (i in 1:mz) {
    for (j in 1:mz) {
      for (t in 1:(T - 1)) {
        myGamma[i, j, t] = alp[t, i] * A[i, j] * B[j, x[t + 1]] * beta[t + 1, j]
      }
    }
  }
  
  A = rowSums(myGamma, dims = 2)
  A = A / rowSums(A)
  
  tmp = apply(myGamma, c(1, 3), sum)  # mz-by-(T-1)
  tmp = cbind(tmp, colSums(myGamma[, , T - 1]))
  
  for (l in 1:mx) {
    B[, l] = rowSums(tmp[, which(x == l)])
  }
  
  B = B / rowSums(B)
  
  return(list(A = A, B = B))
}

myBW = function (x, A, B, w, n.iter = 100) {
  # Input:
  # x: T-by-1 observation sequence
  # A: initial estimate for mz-by-mz transition matrix
  # B: initial estimate for mz-by-mx emission matrix
  # w: initial estimate for mz-by-1 initial distribution over Z_1
  # Output MLE of A and B; we do not update w
  # list(A = A, B = B, w = w)
  
  for (i in 1:n.iter) {
    update.para = BW.onestep(x, A, B, w)
    A = update.para$A
    B = update.para$B
  }
  
  return (list(A = A, B = B, w = w))
}

myViterbi = function (x, A, B, w) {
  states = rownames(A) # A and B
  T = length(x) # Observation in the orginal code
  mz = length(states)
  v = array(NA, c(mz, T))
  
  dimnames(v) = list(states = states, index = 1:T)
  
  # Initialize the values of v to their initial values given the
  # initial distribution of Z and the emission matrix
  # last row of delta
  for (state in states)
  {
    v[state, 1] = log(w[state] * B[state, x[1]])
  }
  
  # Determine the next state of v given the
  # last state and the emission matrix
  for (k in 2:T)
  {
    for (state in states)
    {
      maxi = NULL
      for (previousState in states)
      {
        temp = v[previousState, k - 1] + log(A[previousState, state])
        maxi = max(maxi, temp)
      }
      v[state, k] = log(B[state, x[k]]) + maxi
    }
  }
  
  viterbiPath = rep(NA, T)
  
  for (state in states)
  {
    if (max(v[, T]) == v[state, T])
    {
      viterbiPath[T] = state
      break
    }
  }
  
  # Determine the most likely sequence of states
  for (k in (T - 1):1)
  {
    for (state in states)
    {
      if (max(v[, k] + log(A[, viterbiPath[k + 1]])) == 
          v[state, k] + log(A[state, viterbiPath[k + 1]]))
      {
        viterbiPath[k] = state
        break
      }
    }
  }
  
  
  return(viterbiPath)
}

data = read.csv("Coding3_HMM_Data.csv")

mz = 2
mx = 3
states = levels(data$Z)
symbols = as.factor(data$X)
ini.A = matrix(1, mz, mz)
ini.A = ini.A / rowSums(ini.A)
ini.B = matrix(1:6, mz, mx)
ini.B = ini.B / rowSums(ini.B)
ini.w = c(1 / 2, 1 / 2)

myout = myBW(data$X, ini.A, ini.B, ini.w, n.iter = 100)

rownames(myout$A) = levels(data$Z)
colnames(myout$A) = levels(data$Z)
rownames(myout$B) = levels(data$Z)
names(myout$w) = levels(data$Z)
myout.Z = myViterbi(data$X, myout$A, myout$B, myout$w)

write.table(myout.Z,
            file = "Coding3_HMM_Viterbi_Output.txt",
            row.names = FALSE,
            col.names = FALSE)