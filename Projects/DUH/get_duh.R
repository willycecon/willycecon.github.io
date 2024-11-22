get_duh <- function(data, groups, evepar = 2, round = 4){
require(tidyverse)
  
    if (evepar <= 1){
    stop("Error: Evepar must be strictly greater than 1 for convexity")
  }
  if (length(groups) == 1){
    names <- variable.names(data)
    g <- str_detect(names,
                    as.character(groups))
    G <- sum(g)
    g2 <- names[g]
    subdata <- data[,g] 
    
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
      mutate(ptilde = p/sum(p),
             psi = (ptilde - 1/(G-1))^evepar,
             p2 = p^2,
             plp = (-1) * p * log(p),
             plp = ifelse(is.nan(plp), 0, plp)) %>%
      group_by(id, rmajor) %>%
      summarise(psi = 1 - sum(psi)^(1/evepar),
                HHI = sum(p2),
                SE = sum(plp)) %>%
      mutate(DUH = -(log(rmajor)/log(G))*psi,
             DUH = ifelse(is.nan(DUH), 0, DUH),
             HHI = HHI + rmajor^2,
             GSI = 1 - HHI,
             SE = (-1)*rmajor*log(rmajor) + SE) %>%
      select(id, DUH, GSI, SE,rmajor)
    
    data <- bind_cols(data,
                      rmajor = calc$rmajor,
                      DUH = calc$DUH,
                      GSI = calc$GSI,
                      SE  = calc$SE)  
    print(paste("Universe has", G, "categories",sep=" "))
    print("Version R10b")
    print(paste("Using ",evepar, "-metric", collapse="", sep = ""))
    return(data)
    
  }
}