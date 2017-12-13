template<class I, class O, class P>
O& parameter_advising(I* instance, vector<P> parameters){
  cout << "PA step 0 of " << parameters.size() << endl;
  O* best = instance->solve(parameters[0]);
  cout << "Done" << endl;
  if(parameters[0].smooth) best->improve();
  cout << "Pre-Loop" << endl;
  for(int i=1; i<parameters.size(); i++){
    cout << "PA step " << i << " of " << parameters.size() << endl;
    O* current = instance->solve(parameters[i]);
    if(parameters[i].smooth) current->improve();
    if((*current) > (*best)){
      O* old_best = best;
      best = current;
      delete old_best;
    }else{
      delete current;
    }
  }
  cout << "Post-Loop" << endl;
  return (*best);
}
