tree 	= rncl::read_newick_phylo(file = "tree/full_length.tre")
traits = read.csv("datasets/full_length_reads.csv")
Ntips = length(tree$tip.label)

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

methods = c("MK","MPR","WSCP","SA","EP","PIC")
performance = data.frame()
for (j in 1:6){
  method = methods[j]
  for (i in 0:4){
    labels = data.frame(tree$tip.label)
    colnames(labels) = "accession"
    test = traits[(3904*i+1):(3904*(i+1)),c(2,5)]
    train = traits[,c(2,5)]
    train$copy_number[(3904*i+1):(3904*(i+1))] = NA
    labels = dplyr::left_join(labels,train,by = "accession")
    estimated_states = phylo_estimate(tree = tree,tip_states = labels$copy_number,
                       method = method,Ntips = Ntips,
                       Nstates = 21)
    names(estimated_states) = labels$accession
    predicted = estimated_states[names(estimated_states)%in%test$accession]
    predicted = data.frame(predicted)
    predicted$accession = row.names(predicted)
    test = dplyr::left_join(test,predicted,by = "accession")
    rmse = sqrt(mean((test$copy_number-test$predicted)^2))
    performance[i+1,j]=rmse
  }
}
colnames(performance) = methods

write.csv(performance,"performance/other_algo_full_length.csv",row.names = F)


