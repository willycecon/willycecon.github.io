get_duh <- function(data, groups, evepar = 2, round = 4, show.groups = TRUE, prefix = NULL){
  require(tidyverse)
  
  # data should be a dataframe
  # groups should be a length-1 string that identifies the prefixes of all variables of interests.
  # evepar is the eveness parameter that needs to be strictly greater than 1.
  # round is the rounding parameter for how many digits after 0.
  # show.groups = TRUE means the group names will be shown at the end.
  # prefix takes a length 1 string that adds to the varialbe names of DUH, GSI, and SE.
  if (evepar <= 1){
    stop("Error: evepar must be strictly greater than 1 for convexity")
  }
  if (length(groups) > 1){
    g2 <- groups
    G <- length(groups)
    subdata <- data[,groups]
  }else{
    names <- variable.names(data)
    g <- str_detect(names, regex(paste("^",groups,sep="")))
    G <- sum(g)
    g2 <- names[g]
    subdata <- data[,g] 
  }  
  #Double check groupings are correct
  check <- subdata %>%
    rowwise() %>%
    mutate(rsum = sum(!!!rlang::syms(g2)) %>%
             round(0)-1) %>%
    select(rsum) %>%
    as.matrix() %>%
    as.vector() %>%
    sum()
  if(check != 0){
    # stop("Total Group Proportions Do Not Add Up To 1")
    print("Total Group Proportions Do Not Add Up To 1, automatic scaling applied")
    subdata <- subdata %>%
      rowwise() %>%
      mutate(rsum = sum(!!!rlang::syms(g2))) %>%
      mutate(across(everything(), ~ . / rsum)) %>%
      select(-rsum)
    
    
  }else{
    print("Check: Group Proportions Add Up to 1")
    
    subdata <- data[,g]
  }  
  factor <- log(evepar)
  calc <- subdata %>%
    ungroup() %>%
    mutate(id = 1:n()) %>% 
    pivot_longer(cols = starts_with(as.character(groups)),
                 names_to = "group",
                 values_to = "p") %>%
    group_by(id) %>%
    arrange(id,desc(p)) %>%
    mutate(index = 1:n(),
           rmajor = max(p)) %>%
    filter(index > 1) %>%
    mutate(ptilde = factor*(p/sum(p)),
           psi = abs(ptilde - factor/(G-1))^evepar,
           p2 = p^2,
           plp = (-1) * p * log(p),
           plp = ifelse(is.nan(plp), 0, plp)) %>%
    group_by(id, rmajor) %>%
    summarise(psi = 1 - (sum(psi)^(1/evepar))/(factor),
              HHI = sum(p2),
              SE = sum(plp)) %>%
    mutate(DUH = -(log(rmajor)/log(G))*psi,
           DUH = ifelse(is.nan(DUH), 0, DUH),
           evenness = psi,
           HHI = HHI + rmajor^2,
           GSI = 1 - HHI,
           SE = (-1)*rmajor*log(rmajor) + SE) %>%
    select(id, DUH, GSI, SE, rmajor, evenness)
  data <- bind_cols(data,
                    rmajor = calc$rmajor,
                    evenness = calc$evenness,
                    DUH = calc$DUH,
                    GSI = calc$GSI,
                    SE  = calc$SE)
  
  if (!is.null(prefix)){
    name <- names(data)
    name_length <- length(name)
    name[(name_length-2):name_length] <- c(paste(prefix,"DUH",sep=""),
                                           paste(prefix,"GSI",sep=""),
                                           paste(prefix,"SE",sep=""))
    names(data) <- name
  }
  
  print(paste("Universe has", G, "categories",sep=" "))
  if(show.groups == TRUE){
    print(paste(g2, collapse = ", "))
  }
  # print("Version R10b")
  print(paste("Using ",evepar, "-metric", collapse="", sep = ""))
  return(data)
}


