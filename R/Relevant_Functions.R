# set up group information
# m: group multiplier, default (missing) is the square root of group size of remaining features
setupG <- function(group, m){
  group.factor <- factor(group)
  if (any(levels(group.factor) == '0')) {
    g <- as.integer(group.factor) - 1
    lev <- levels(group.factor)[levels(group.factor) != '0']
  } else {
    g <- as.integer(group.factor)
    lev <- levels(group.factor)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (missing(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    TRY <- try(as.integer(group) == g)
    if (inherits(TRY, 'try-error') || any(!TRY)) stop('Attempting to set group.multiplier is ambiguous if group is not a factor', call. = FALSE)
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups", call. = FALSE)
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
    if (any(m < 0)) stop('group.multiplier cannot be negative', call.=FALSE)
  }
  # "g" contains the group index of each column, but convert "character" group name into integer
  structure(g, levels = lev, m = m)
}

# remove constant columns if necessary
subsetG <- function(g, nz) { # nz: index of non-constant features
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  new <- g[nz] # only include non-constant columns
  dropped <- setdiff(g, new)  # If the entire group has been dropped
  if (length(dropped) > 0) {
    lev <- lev[-dropped] # remaining group
    m <- m[-dropped]
    group.factor <- factor(new) #remaining group factor
    new <- as.integer(group.factor) - 1 * any(levels(group.factor) == '0')  #new group index
  }
  structure(new, levels = lev, m = m)
}

# reorder group index of features
reorderG <- function(g, m) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g == 0)) {
    g <- as.integer(relevel(factor(g), "0")) - 1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf) == "0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels = lev, m = m, ord = ord, ord.inv = ord.inv, reorder = reorder)
}

#Feather-level standardization
standardize.Z <- function(Z){
  mysd <- function(z){
    sqrt(sum((z - mean(z))^2)/length(z))
  }
  new.Z <- scale(as.matrix(Z), scale = apply(as.matrix(Z), 2, mysd))
  center.Z <- attributes(new.Z)$`scaled:center`
  scale.Z <- attributes(new.Z)$`scaled:scale`
  new.Z <- new.Z[, , drop = F]
  res <- list(new.Z = new.Z, center.Z = center.Z, scale.Z = scale.Z)
  return(res)
}

## converting standardized betas back to original variables
unstandardize <- function(beta, gamma, std.Z){
  original.beta <- matrix(0, nrow = length(std.Z$scale), ncol = ncol(beta))
  original.beta[std.Z$nz, ] <- beta / std.Z$scale[std.Z$nz]  # modified beta
  original.gamma <- t(apply(gamma, 1, function(x) x - crossprod(std.Z$center, original.beta))) # modified intercepts (gamma)
  return(list(gamma = original.gamma, beta = original.beta))
}


# Group-level orthogonalization (column in new order, from group_0 to group_max)
orthogonalize <- function(Z, group) {
  z.names <- colnames(Z)
  n <- nrow(Z)
  J <- max(group)
  QL <- vector("list", J)
  orthog.Z <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  colnames(orthog.Z) <- z.names
  # unpenalized group will not be orthogonalized
  orthog.Z[, which(group == 0)] <- Z[, which(group == 0)]

  # SVD and generate orthogonalized X
  for (j in seq_along(integer(J))) {
    ind <- which(group == j)
    if (length(ind) == 0) { # skip 0-length group
      next
    }
    SVD <- svd(Z[, ind, drop = FALSE], nu = 0)  # Q matrix (orthonormal matrix of eigenvectors)
    r <- which(SVD$d > 1e-10)  #remove extremely small singular values
    QL[[j]] <- sweep(SVD$v[, r, drop = FALSE], 2, sqrt(n)/SVD$d[r], "*") # Q * Lambda^{-1/2}
    orthog.Z[, ind[r]] <- Z[, ind] %*% QL[[j]]  # group orthogonalized X, where (X^T * X)/n = I
  }
  nz <- !apply(orthog.Z == 0, 2, all)  #find all zero
  orthog.Z <- orthog.Z[, nz, drop = FALSE]
  attr(orthog.Z, "QL") <- QL
  attr(orthog.Z, "group") <- group[nz]
  return(orthog.Z)
}

# convert orthogonalized beta back to original scales
unorthogonalize <- function(beta, Z, group) {
  ind <- !sapply(attr(Z, "QL"), is.null)
  QL <- Matrix::bdiag(attr(Z, "QL")[ind]) #block diagonal matrix
  if (sum(group == 0) > 0){ #some groups are unpenalized
    ind0 <- which(group==0)
    original.beta <- as.matrix(rbind(beta[ind0, , drop = FALSE], QL %*% beta[-ind0, , drop = FALSE]))
  } else {  # all groups are penalized
    original.beta <- as.matrix(QL %*% beta)
  }
  return(original.beta)
}


# standardize + orthogonalize covariate matrix
newZG.Std.grplasso <- function(data, Z.char, g, m){
  Z <- as.matrix(data[, Z.char, drop = F])
  if (any(is.na(Z))){
    stop("Missing data (NA's) detected in covariate matrix!", call. = FALSE)
  }
  if (length(g) != ncol(Z)) {
    stop ("Dimensions of group is not compatible with Z", call. = FALSE)
  }
  G <- setupG(g, m) # setup group
  std <- standardize.Z(Z)
  std.Z <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  
  small_scales <- which(scale <= 1e-6)
  if (length(small_vals) > 0) {
    stop(
      paste0(
        "The following variables have (near) constant columns: ",
        paste(names(scale)[small_vals], collapse = ", ")
      )
    )
  }
  
  
  nz <- which(scale > 1e-6)   # non-constant columns
  if (length(nz) != ncol(Z)) {
    std.Z <- std.Z[, nz, drop = F]
    G <- subsetG(G, nz)
  }
  # Reorder groups
  G <- reorderG(G, attr(G, 'm'))
  if (attr(G, 'reorder')){
    std.Z <- std.Z[, attr(G, 'ord')]
  }
  # Group-level orthogonalization
  std.Z <- orthogonalize(std.Z, G)
  g <- attr(std.Z, "group")
  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- sqrt(table(g[g != 0]))
  }
  res <- list(std.Z = std.Z, g = g, m = m, reorder = attr(G, 'reorder'), nz = nz,
              ord.inv = attr(G, 'ord.inv'), center = center, scale = scale)
  return(res)
}

# Only orthogonalize covariate matrix
newZG.Unstd.grplasso <- function(data, Z.char, g, m){
  Z <- as.matrix(data[, Z.char, drop = F])
  if (any(is.na(Z))){
    stop("Missing data (NA's) detected in covariate matrix!", call. = FALSE)
  }
  if (length(g) != ncol(Z)) {
    stop ("Dimensions of group is not compatible with Z", call. = FALSE)
  }
  G <- setupG(g, m)
  mysd <- function(x){
    sqrt(sum((x - mean(x))^2)/length(x))
  }
  scale <- apply(as.matrix(Z), 2, mysd)
  
  
  small_scales <- which(scale <= 1e-6)
  if (length(small_vals) > 0) {
    stop(
      paste0(
        "The following variables have (near) constant columns: ",
        paste(names(scale)[small_vals], collapse = ", ")
      )
    )
  }
  
  nz <- which(scale > 1e-6) #remove constant columns
  if (length(nz) != ncol(Z)) {
    std.Z <- Z[, nz, drop = F]
    G <- subsetG(G, nz)
  } else {
    std.Z <- Z
  }
  G <- reorderG(G, attr(G, 'm'))
  if (attr(G, 'reorder')){
    std.Z <- std.Z[, attr(G, 'ord')]
  }
  std.Z <- orthogonalize(std.Z, G)
  g <- attr(std.Z, "group")
  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- sqrt(table(g[g != 0]))
  }
  res <- list(std.Z = std.Z, g = g, m = m, reorder = attr(G, 'reorder'),
              ord.inv = attr(G, 'ord.inv'), nz = nz)
  return(res)
}


## handle multiple columns of Y
newY <- function(data, Y.char){
  y <- data[, Y.char]
  if (is.data.frame(y)){
    y <- as.matrix(y)
  }
  if (is.matrix(y)) {
    d <- dim(y)
    y <- t(y)
  } else {
    d <- c(length(y), 1)
  }

  # Convert fuzzy binomial data
  if (typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) stop("Outcome is not binary", call.=FALSE)
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(Y = ", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])  #convert to 0 & 1
      if (d[2] > 1) {
        attr(y, "dim") <- d
      }
    }
  }

  # Convert to double, if necessary
  if (typeof(y) != "double") {
    tryCatch(storage.mode(y) <- "double", warning = function(w) {stop("Y must be numeric or able to be coerced to numeric", call. = FALSE)})
  }
  if (any(is.na(y))){
    stop("Missing data (NA's) detected in outcome Y!", call. = FALSE)
  }

  # Handle multi
  if (is.matrix(y)) {
    if (ncol(y) > 1) {
      if (is.null(colnames(y))){
        paste("Y", 1:ncol(y), sep = "")
      }
    }
    attributes(y) <- NULL
  }

  attr(y, "m") <- d[2]
  return(y)
}

