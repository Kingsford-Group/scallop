Template<class I, class O, class P>
O parameter_advising(I instance, vector<P> parameters){
  O best = I.solve(parameters[0]);
  if(P.smooth) best.improve();
  for(int i=1; i<parameters.size(); i++){
    O current = I.solve(parameters[i]);
    if(P.smooth) current.improve();
    if(current > best) best = current;
  }
  return current;
}
