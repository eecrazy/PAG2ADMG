library(pcalg)
library(ggm)

#   Function that makes a random pag given a specific number of nodes, latent variable indices
#   probability of edge between any pair of nodes and an alpha value for the independence test
make_rdm_pag <- function(num_nodes, latent_idx, prob, alpha){
    g <- randomDAG(n = num_nodes, prob = prob)
    true.corr <- cov2cor(trueCov(g))
    true.pag <- dag2pag(suffStat = list(C= true.corr, n= 10^9),
                        indepTest= gaussCItest, graph=g, L=latent_idx, 
                        alpha= alpha)
    return(true.pag)
}

#   Function that takes an adjacency matrix of a dag and latent indices and returns a pag
make_pag <- function(dag_amat, latent_idx){
    dag <- new('graphAM', adjMat = dag_amat, edgemode = 'directed')
    dag <- as(dag, 'graphNEL')
    true.corr <- cov2cor(trueCov(dag))
    true.pag <- dag2pag(suffStat = list(C= true.corr, n= 10^9),
                        indepTest= gaussCItest, graph=dag, L=latent_idx, 
                        alpha= 0.05)
    return(true.pag)
}

#   Function that takes in a pag and enumerates all possible ADMGs on those nodes based on the
#   intelligent enumeration method provided in eval_possibilities_ggm
pag2admg_verbose <- function(pag){
    pag_amat <- pag@amat
    admg_list <- list()
    poss_list <- vector() #nothing, undirected, directed, 
    node_combos <- combn(rownames(pag_amat), 2)
    node_combos_names <- vector()
    for(i in 1:ncol(node_combos)){
        node_names <- node_combos[,i]
        poss_list <- append(poss_list, list(eval_possibilities_ggm(pag_amat, node_names)))
        node_combos_names <- c(node_combos_names, paste(node_names[1], node_names[2], sep = '*'))
    }
    
    poss_grid <- expand.grid(poss_list)
    colnames(poss_grid) <- node_combos_names
    
    for(i in 1:nrow(poss_grid)){
        new_amat <- build_valid_admg_amat(poss_grid[i,], ncol(pag_amat), pag_amat)
        if(is.null(new_amat) == FALSE){
            admg_list <- append(admg_list, list(new_amat))
        }
    }
    return(admg_list)
}

#   Function that takes in a pag and enumerates all possible ADMGs on those nodes based on the
#   brute force enumeration method provided in eval_possibilities_ggm_brute
pag2admg_verbose_brute <- function(pag){
    pag_amat <- pag@amat
    admg_list <- list()
    poss_list <- vector() #nothing, undirected, directed, 
    node_combos <- combn(rownames(pag_amat), 2)
    node_combos_names <- vector()
    for(i in 1:ncol(node_combos)){
        node_names <- node_combos[,i]
        poss_list <- append(poss_list, list(eval_possibilities_ggm_brute(pag_amat, node_names)))
        node_combos_names <- c(node_combos_names, paste(node_names[1], node_names[2], sep = '*'))
    }
    
    poss_grid <- expand.grid(poss_list)
    colnames(poss_grid) <- node_combos_names
    
    for(i in 1:nrow(poss_grid)){
        new_amat <- build_valid_admg_amat(poss_grid[i,], ncol(pag_amat), pag_amat)
        #print(new_amat)
        if(is.null(new_amat) == FALSE){
            admg_list <- append(admg_list, list(new_amat))
        }
    }
    return(admg_list)
}

#   Function that takes an admg adjacency matrix and converts it into a DAG with latent variable confounders
#   in every location that has a bi-directed edge, then converts it to a PAG and returns that PAG
admg2dag2pag <- function(admg_amat){
    num_nodes <- ncol(admg_amat)
    for (i in 1:num_nodes){
        bidirected_idx <- which(admg_amat[,i] == 100 | admg_amat[,i] == 101)
        num_bidirected <- length(bidirected_idx)
        if(num_bidirected > 0){
            for (j in 1:length(bidirected_idx)){
                node_two <- bidirected_idx[j]
                #print(admg_amat[node_two,i])
                admg_amat[node_two,i] <- admg_amat[node_two,i] - 100
                admg_amat[i,node_two] <- admg_amat[i,node_two] - 100
                new_node_row <- rep(0, ncol(admg_amat))
                new_node_row[i] <- 1
                new_node_row[node_two] <- 1
                new_node_col <- rep(0, ncol(admg_amat) + 1)
                admg_amat <- rbind(admg_amat, new_node_row)
                admg_amat <- cbind(admg_amat, new_node_col)
                colnames(admg_amat) <- rownames(admg_amat) <- 1:nrow(admg_amat)
            }
        }
    }
    
    if(ncol(admg_amat) > num_nodes){
        latent_idx <- (num_nodes+1):ncol(admg_amat)
    } else {
        latent_idx <- NULL
    }
    
    new_dag_amat <- admg_amat
    #print(new_dag_amat)
    dag <- new('graphAM', adjMat = new_dag_amat, edgemode = 'directed')
    dag <- as(dag, 'graphNEL')
    #print(dag)
    true.corr <- cov2cor(trueCov(dag))
    true.pag <- dag2pag(suffStat = list(C= true.corr, n= 10^9),
                        indepTest= gaussCItest, graph=dag, L=latent_idx, 
                        alpha= 0.9999)
    return(true.pag)
}

#   Function that takes a list of possible edges, number of nodes, and the pag amat
#   Returns a valid admg if one exists based on that information
build_valid_admg_amat <- function(edge_info, num_nodes, pag_amat){
    amat <- matrix(rep(0,num_nodes^2), nrow = num_nodes, ncol = num_nodes)
    rownames(amat) <- rownames(pag_amat)
    colnames(amat) <- colnames(pag_amat)
    column_names <- colnames(edge_info)
    for(i in 1:length(column_names)){
        col_name <- column_names[i]
        star_idx <- which(strsplit(col_name, '')[[1]] == '*')
        node_one <- substr(col_name, 1, star_idx-1)
        node_two <- substr(col_name, star_idx+1, nchar(col_name))
        amat[node_one, node_two] <- edge_info[[i]][[1]][1]
        amat[node_two, node_one] <- edge_info[[i]][[1]][2]
    }
    
    if(isADMG(amat)){
        return(amat)
    } else {
        return(NULL)
    }
    
}

#   Function that enumerates the nodes based on the pcalg definition of edges (Not used)
eval_possibilities_pcalg <- function(amat, node_names){
    node_one <- node_names[1]
    node_two <- node_names[2]
    sig_two <- amat[node_one, node_two]
    sig_one <- amat[node_two, node_one]
    
    #   Check if 0 -- 0
    if(sig_one == 1 & sig_two == 1){
        return(list(c(2,3), c(3,2)))
    }
    
    if(sig_one == 1 & sig_two == 2){
        return(list(c(0,0), c(2,2), c(3,2)))
    }
    
    if(sig_one == 2 & sig_one == 1){
        return(list(c(0,0), c(2,2), c(2,3)))
    }
    
    if(sig_one == 2 & sig_two == 2){
        return(list(c(2,2)))
    }
    
    if(sig_one == 0 & sig_two == 0){
        return(list(c(0,0)))
    }
    
    if(sig_one == 3 & sig_two == 2){
        return(list(c(3,2)))
    }
    
    if(sig_one == 2 & sig_two == 3){
        return(list(c(2,3)))
    }
}

#   Function that enumerates the nodes based on the ggm definition of edges and does this in the
#   intelligent manner provided in the report
eval_possibilities_ggm <- function(amat, node_names){
    node_one <- node_names[1]
    node_two <- node_names[2]
    sig_two <- amat[node_one, node_two]
    sig_one <- amat[node_two, node_one]
    
    #   Check if 0 -- 0
    if(sig_one == 1 & sig_two == 1){
        return(list(c(1,0), c(0,1), c(0,0), c(100,100), c(101,100), c(100,101)))
    }
    
    if(sig_one == 1 & sig_two == 2){
        return(list(c(1,0), c(0,0), c(100,100), c(101,100)))
    }
    
    if(sig_one == 2 & sig_two == 1){
        return(list(c(0,1), c(0,0), c(100,100), c(100,101)))
    }
    
    if(sig_one == 2 & sig_two == 2){
        return(list(c(0,0), c(100,100)))
    }
    
    if(sig_one == 0 & sig_two == 0){
        return(list(c(0,0)))
        #return(list(c(1,0), c(0,1), c(0,0), c(100,100)))
    }
    
    if(sig_one == 3 & sig_two == 2){
        return(list(c(1,0)))
    }
    
    if(sig_one == 2 & sig_two == 3){
        return(list(c(0,1)))
    }
}

#   Function that enumerates the nodes based on the ggm definition of edges and does this in the
#   brute force manner
eval_possibilities_ggm_brute <- function(amat, node_names){
    node_one <- node_names[1]
    node_two <- node_names[2]
    sig_two <- amat[node_one, node_two]
    sig_one <- amat[node_two, node_one]
    
    #   Check if 0 -- 0
    if(sig_one == 1 & sig_two == 1){
        return(list(c(1,0), c(0,1), c(0,0), c(100,100), c(101,100), c(100,101)))
    }
    
    if(sig_one == 1 & sig_two == 2){
        return(list(c(1,0), c(0,1), c(0,0), c(100,100), c(101,100), c(100,101)))
    }
    
    if(sig_one == 2 & sig_two == 1){
        return(list(c(1,0), c(0,1), c(0,0), c(100,100), c(101,100), c(100,101)))
    }
    
    if(sig_one == 2 & sig_two == 2){
        return(list(c(1,0), c(0,1), c(0,0), c(100,100), c(101,100), c(100,101)))
    }
    
    if(sig_one == 0 & sig_two == 0){
        return(list(c(0,0)))
        #return(list(c(1,0), c(0,1), c(0,0), c(100,100)))
    }
    
    if(sig_one == 3 & sig_two == 2){
        return(list(c(1,0), c(0,1), c(0,0), c(100,100), c(101,100), c(100,101)))
    }
    
    if(sig_one == 2 & sig_two == 3){
        return(list(c(1,0), c(0,1), c(0,0), c(100,100), c(101,100), c(100,101)))
    }
}

#   Function that takes in a pag and returns the admg list that describes this using the brute
#   force heuristic
brute_force_function <- function(pag){
    admg_list <- pag2admg_verbose_brute(pag)
    parsed_admg_list <- list()
    for(i in 1:length(admg_list)){
        temp_admg <- admg_list[[i]]
        temp_pag <- admg2dag2pag(temp_admg)
        if(all(temp_pag@amat == pag@amat)){
            parsed_admg_list <- append(parsed_admg_list, list(temp_admg))
        }
    }
    return(parsed_admg_list)
}

#   Function that takes in a pag and returns the admg list that describes this using the intelligent
#   heuristic provided in the paper
pag2admg <- function(pag){
    admg_list <- pag2admg_verbose(pag)
    print(length(admg_list))
    parsed_admg_list <- list()
    for(i in 1:length(admg_list)){
        temp_admg <- admg_list[[i]]
        temp_pag <- admg2dag2pag(temp_admg)
        if(all(temp_pag@amat == pag@amat)){ #PAG Equivalence Check
            parsed_admg_list <- append(parsed_admg_list, list(temp_admg))
        }
    }
    return(parsed_admg_list)
}

#   Function that allows enumeration of all edge types based on a pag and a list of admgs
enumerate_all_resulting_edges <- function(pag_amat, admg_amat_list){
    enumeration_amat <- matrix(nrow = nrow(pag_amat), ncol = ncol(pag_amat))
    for(i in 1:nrow(pag_amat)){
        for(j in 1:ncol(pag_amat)){
            choices <- vector()
            for(k in 1:length(admg_amat_list)){
                val <- admg_amat_list[[k]][i,j]
                if(val %in% choices == FALSE){
                    choices <- c(choices, val)
                }
            }
            enumeration_amat[i,j] <- paste(choices, collapse = ' ')
        }
    }
    return(enumeration_amat)
}

#   Function that plots all admgs in the external plot window
plot_all_admgs <- function(admg_list){
    for(i in 1:length(admg_list)){
        plotGraph(admg_list[[i]], layout = layout.circle)
    }
}

#   Function that plots all admgs in R using pcalgs definition
plot_all_admgs_pcalg <- function(admg_list){
    for(i in 1:length(admg_list)){
        plot(convert_admg_ggm_to_pcalg_to_plot(admg_list[[i]]))
    }
}

#   Takes a pag adjacency matrix and creates a pag object
make_pag_from_amat <- function(pag_amat){
    rownames(pag_amat) <- colnames(pag_amat) <- 1:nrow(pag_amat)
    pag <- new('fciAlgo', amat = pag_amat)
    #pag <- as(pag, 'graphNEL')
    return(pag)
}

#   Converts an admg_ggm_amat into an amat graphable in pcalg in R
convert_admg_ggm_to_pcalg_to_plot <- function(admg_ggm_amat){
    temp <- admg_ggm_amat + 10
    for(i in 1:nrow(temp)){
        bi_idx <- which(temp[i,] == 110 || temp[i,] == 111)
        dir_idx <- which(temp[i,] == 11 || temp[i,] == 111)
        admg_ggm_amat[i, bi_idx] <- 2
        admg_ggm_amat[i, dir_idx] <- 2
        admg_ggm_amat[dir_idx, i] <- 3
    }
    pag_amat <- admg_ggm_amat
    rownames(pag_amat) <- colnames(pag_amat) <- 1:nrow(pag_amat)
    pag <- new('fciAlgo', amat = pag_amat)
    return(pag)
}

#   Checks if the input adj matrix refers to a valid ADMG
isADMG <- function(amat){
    ### check is if a graph is an ADMG
    comp <- unmakeMG(amat)
    ug <- comp$ug; dag <- comp$dg; bg <- comp$bg  
    out <- TRUE
    # if(any(amat > 100)){  
    #     warning("There are double edges.")
    #     out <- FALSE
    # }
    if(!isAcyclic(dag)){
        warning("Not acyclic.")
        out <- FALSE
    }
    out 
}  
