template <typename T>
Heap<T>::Heap() {
  array.clear();
}

template <typename T>
void Heap<T>::insert(T item) {
  array.push_back(item);
  BubbleUp(array.size()-1);
}

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



template<typename T>
void Heap<T>::print() {

  for(unsigned i=0; i<array.size(); ++i) {
    std::cout << array[i] << std::endl;
  }
}
