# Align latent positions across time to prevent arbitrary sign/rotation flips
# Unipartite: U is [n, R, T]
align_over_time_unip <- function(U) {
  if (is.null(U) || length(dim(U)) != 3L) return(U)
  T <- dim(U)[3]
  if (T < 2L) return(U)
  U_al <- U
  for (t in 2:T) {
    M <- crossprod(U_al[,,t], U_al[,,t-1])  # R x R
    sv <- base::svd(M)
    Ropt <- sv$u %*% t(sv$v)
    U_al[,,t] <- U_al[,,t] %*% Ropt
  }
  U_al
}

# Bipartite: U [nA, RA, T], V [nB, RB, T], G [RA, RB]
align_over_time_bip <- function(U, V, G) {
  if (is.null(U) || is.null(V) || length(dim(U)) != 3L || length(dim(V)) != 3L) {
    return(list(U = U, V = V, G = G))
  }
  T <- dim(U)[3]
  if (T < 2L) return(list(U = U, V = V, G = G))
  U_al <- U; V_al <- V; G_al <- G
  for (t in 2:T) {
    MU <- crossprod(U_al[,,t], U_al[,,t-1])
    MV <- crossprod(V_al[,,t], V_al[,,t-1])
    su <- svd(MU); sv <- svd(MV)
    RU <- su$u %*% t(su$v)  # RA x RA
    RV <- sv$u %*% t(sv$v)  # RB x RB
    U_al[,,t] <- U_al[,,t] %*% RU
    V_al[,,t] <- V_al[,,t] %*% RV
    # keep U G V' invariant over rotations
    G_al <- t(RU) %*% G_al %*% RV
  }
  list(U = U_al, V = V_al, G = G_al)
}

# Estimate AR(1)-like persistence directly from aligned latent series
estimate_rho_from_U <- function(U) {
  if (is.null(U) || length(dim(U)) != 3L) return(NA_real_)
  T <- dim(U)[3]
  if (T < 2L) return(NA_real_)
  num <- 0; den <- 0
  for (t in 2:T) {
    Ut  <- U[,,t]
    Ut1 <- U[,,t-1]
    num <- num + sum(Ut * Ut1, na.rm = TRUE)
    den <- den + sum(Ut1 * Ut1, na.rm = TRUE)
  }
  if (den > 0) num / den else NA_real_
}