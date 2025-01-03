## 4. Dynamic Data Structure

1. When dynamic analysis is performed, the results from static analysis need to be transfered to a data container so that it is not destroyed at the end of a frame. The contrainer is a list which has a length same as the frame interval. At the beginning of dynamic analysis run, the container will be initialized. For each frame, elements in the list are pushed from last to first, so that the last position is empty and can be used to store the latest data. This function is named `collect_data()`.

2. After the data is stored in the contrainer, the function named `compare_data()` can be called to compare the data in the container at different frames. In summary, the `collect_data()` should have access to the static analysis modules, and the `compare_data()` should have access to the dynamic data container.

3. A draw back is that, since the static analysis results are destroyed (by deallocate the variable), the value is transferred to dynamic data container by value. This means the speed can be slow when the analysis results are huge. A way to optimize it is to use pointer to point to the results of static analysis, nullify the pointer and let the dynamic data container to point to the results. In this way the data will be transferred by address.


