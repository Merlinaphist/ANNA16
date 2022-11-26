phylo_estimate = function(tree,tip_states,method,Ntips,Nstates){
  if (method == "MK"){
    hsp = castor::hsp_mk_model(tree,tip_states,Nstates)
  }else if(method == "MPR"){
    hsp = castor::hsp_max_parsimony(tree,tip_states,Nstates)
  }else if(method == "WSCP"){
    hsp = castor::hsp_squared_change_parsimony(tree,tip_states,Nstates,weighted=T)
  }else if(method == "SA"){
    hsp = castor::hsp_subtree_averaging(tree,tip_states,Nstates)
  }else if(method == "EP"){
    hsp = castor::hsp_empirical_probabilities(tree,tip_states,Nstates)
  }else if(method == "PIC"){
    hsp = castor::hsp_independent_contrasts(tree,tip_states,Nstates)
  }
  if(method %in% c("WSCP","PIC","SA")){
    estimated_state_values = hsp$states[1:Ntips] # round estimated numeric values to nearest integer (since fractional trait values don't make sense for 16S GCNs)
  }else{
    # find maximum likelihood states based on posterior probabilities
    estimated_state_values = max.col(hsp$likelihoods[1:Ntips,]) # map max-likelihood states back to numeric value
  }
  return(estimated_state_values)
}

load_tree = function(filename){
  
  tree = rncl::read_newick_phylo(file = filename)
  
  TREE_EDGE_LENGTH_EPSILON=NULL
  
  if(is.null(tree$node.label)){
    cat(sprintf("Adding node labels to full tree..\n"))
    tree$node.label = paste("node.", 1:Nnodes, sep = "") # don't use underscores, because some tree readers (e.g. rncl) interpret them as spaces
  }
  
  # avoid zero-length edges
  if(any(tree$edge.length==0)){
    if(is.null(TREE_EDGE_LENGTH_EPSILON)) TREE_EDGE_LENGTH_EPSILON = 0.1*min(tree$edge.length[tree$edge.length>0])
    cat(sprintf("Note: Some edges have length zero, which may break some of the HSP routines. Replacing zero-lengths with a tiny positive length (%g)..\n",TREE_EDGE_LENGTH_EPSILON))
    tree$edge.length[tree$edge.length==0] = TREE_EDGE_LENGTH_EPSILON
  }
  return(tree)
}


traits = read.csv("datasets/full_length_reads.csv")
regions = c("V1-V2","V1-V3","V3-V4","V4","V4-V5","V6-V8","V7-V9")
methods = c("MPR","WSCP","SA","EP","PIC")

all_performance = data.frame()
for (k in 1:5){
  method = methods[k]
  performance = data.frame()
  for (j in 1:7){
    region = regions[j]
    tree = load_tree(paste0("tree/",region,".tre"))
    Ntips = length(tree$tip.label)
    labels = data.frame(tree$tip.label)
    colnames(labels) = "accession"
    tmp_traits = traits[(traits$accession %in% labels$accession),]
    multiplicand = as.integer(nrow(tmp_traits)*0.2)+1
    for (i in 0:4){
      labels = data.frame(tree$tip.label)
      colnames(labels) = "accession"
      if (i == 4){
        test = tmp_traits[(multiplicand*i+1):nrow(tmp_traits),c(2,5)]
        train = tmp_traits[,c(2,5)]
        train$copy_number[(multiplicand*i+1):nrow(tmp_traits)] = NA
      }else{
        test = tmp_traits[(multiplicand*i+1):(multiplicand*(i+1)),c(2,5)]
        train = tmp_traits[,c(2,5)]
        train$copy_number[(multiplicand*i+1):(multiplicand*(i+1))] = NA
      }
      labels = dplyr::left_join(labels,train,by = "accession")
      estimated_states = phylo_estimate(tree,labels$copy_number,method=method,Ntips,21)
      names(estimated_states) = labels$accession
      predicted = estimated_states[names(estimated_states)%in%test$accession]
      predicted = data.frame(predicted)
      predicted$accession = row.names(predicted)
      test = dplyr::left_join(test,predicted,by = "accession")
      rmse = sqrt(mean((test$copy_number-test$predicted)^2))
      performance[i+1,j]=rmse
    }
  }
  colnames(performance) = regions
  performance$model = method
  all_performance[(5*k-4):(5*k),1:8] = performance
}

write.csv(all_performance,"performance/other_algo_subregions.csv",row.names = F)
