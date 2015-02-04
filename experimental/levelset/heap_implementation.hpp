//=========================================
template <typename T>
Heap<T>::Heap() {
  array.clear();
}

//=========================================
template <typename T>
void Heap<T>::insert(T item) {
  array.push_back(item);
  BubbleUp(array.size()-1);
}

//=========================================
template<typename T>
void Heap<T>::BubbleUp(unsigned node) {
  if(node==0)
    return;

  unsigned p = parent(node);
  if(array[node] < array[p]) {
    T tmp = array[node];
    array[node] = array[p];
    array[p] = tmp;
    BubbleUp(p);
  }
  return;
}

//=========================================
template<typename T>
T Heap<T>::extract() {
  
  T retval = array[0];
  array[0] = array[array.size()-1];
  array.pop_back();
  BubbleDown(0);
  
  return retval;
}

//=========================================
template<typename T>
void Heap<T>::BubbleDown(unsigned node) {
  
  unsigned Li, Ri, min;
  Li = 2*node + 1;
  Ri = 2*node + 2;
  min=node;
  
  if(Li < array.size() && array[Li]<array[min])
    min = Li;
  if(Ri < array.size() && array[Ri] < array[min])
    min = Ri;
  
  if(min != node) {
    T tmp = array[node];
    array[node] = array[min];
    array[min] = tmp;
    
    BubbleDown(min);
  }
}

//=========================================
template<typename T>
void Heap<T>::print() {

  for(unsigned i=0; i<array.size(); ++i) {
    std::cout << array[i] << std::endl;
  }
}
