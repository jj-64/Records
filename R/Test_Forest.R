
Forest = function(x, sig=0.05){

## Tree: dataframe of nrow=simulations, returning the decision of each tree[column] at each simulation[row]
tree = as.data.frame(matrix(0,ncol=24,nrow=1))
colnames(tree) = c("CLYD","CLDY","CYLD","CYDL","CDYL","CDLY",
                   "LCYD","LCDY","LYCD","LYDC","LDYC","LDCY",
                   "YCLD","YCDL","YLCD","YLDC","YDLC","YDCL",
                   "DCLY","DCYL","DYCL","DYLC","DLYC","DLCY"
                      )

## Result for the average: what the whole sequence of trees is voting for (simulation[row] and model[col])
result=data.frame(matrix(0,nrow=1,ncol=5))
colnames(result)=c("None","Classical","DTRW","LDM","YN")

############################ Perform simulation ################################

  tree[1,1] = DT_CLYD(X=x, p =sig)
  tree[1,2] = DT_CLDY(X=x, p =sig)
  tree[1,3] = DT_CYLD(X=x, p =sig)
  tree[1,4] = DT_CYDL(X=x, p =sig)
  tree[1,5] = DT_CDYL(X=x, p =sig)
  tree[1,6] = DT_CDLY(X=x, p =sig)

  tree[1,7] = DT_LCYD(X=x, p =sig)
  tree[1,8] = DT_LCDY(X=x, p =sig)
  tree[1,9] = DT_LYCD(X=x, p =sig)
  tree[1,10] = DT_LYDC(X=x, p =sig)
  tree[1,11] = DT_LDYC(X=x, p =sig)
  tree[1,12] = DT_LDCY(X=x, p =sig)

  tree[1,13] = DT_YCLD(X=x, p =sig)
  tree[1,14] = DT_YCDL(X=x, p =sig)
  tree[1,15] = DT_YLCD(X=x, p =sig)
  tree[1,16] = DT_YLDC(X=x, p =sig)
  tree[1,17] = DT_YDLC(X=x, p =sig)
  tree[1,18] = DT_YDCL(X=x, p =sig)

  tree[1,19] = DT_DCLY(X=x, p =sig)
  tree[1,20] = DT_DCYL(X=x, p =sig)
  tree[1,21] = DT_DYCL(X=x, p =sig)
  tree[1,22] = DT_DYLC(X=x, p =sig)
  tree[1,23] = DT_DLYC(X=x, p =sig)
  tree[1,24] = DT_DLCY(X=x, p =sig)

## Summary
  ## Results
  result[1,"None"] = sum(tree[1,]==0)
  result[1,"Classical"] = sum(tree[1,]==1)
  result[1,"DTRW"] = sum(tree[1,]==2)
  result[1,"LDM"] = sum(tree[1,]==3)
  result[1,"YN"] = sum(tree[1,]==4)

  return(list("Tree"=tree,"summary" = result))
}
