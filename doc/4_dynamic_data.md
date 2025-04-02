## 4. Dynamic Data Structure

### General information

1. When dynamic analysis is performed, the results from static analysis need to be transfered to a data container so that it is not destroyed at the end of a frame. The contrainer is a list which has a length same as the frame interval. For different computed data, the data container would also be different. At the beginning of dynamic analysis run, the container will be initialized. At each step, data in the container is pushed from last to first, so that the last position is empty and can be used to store the latest data. This function is named `collect_data()`. For different computed data, the corresponding `collect_xxx()` will be called inside the function. So, following parts need to be adjusted when developing a new method of dynamic analysis:
    - Add corresponding data container in the dynamic data module;
    - Add collect_data function in the main module to push the static analysis data into the container;
    - Add compare_data function in the main module to call the actual analysis function;
    - Write the actual analysis part in its corresponding module.

2. After the data is stored in the contrainer, the function named `compare_data()` can be called to compare the data in the container at different frames. For different computed data, corresponding `compare_xxx()` will be called during comparison. For accessbility, the `collect_data()` should have access to the static analysis modules, and the `compare_data()` should have access to the dynamic data container.

3. A draw back is that, since the static analysis results are destroyed (by deallocate the variable), the value is transferred to dynamic data container by value. This means the speed can be slow when the analysis results are huge. A way to optimize it is to use pointer to point to the results of static analysis, nullify the pointer and let the dynamic data container to point to the results. In this way the data will be transferred by address.

### Specific settings

1. For analysis including:
    - Neighbor change analysis;
    - D2min analysis;
    - Polyhedral neighbor change analysis;

    , there are two kinds of reference configurations. The first is that the reference configuration is always the first frame. The other kind is that the reference configuration has a constant interval to current frame.
