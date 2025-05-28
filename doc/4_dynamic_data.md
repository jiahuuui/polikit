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

### Post process

1. The analysis results will be all included in one log file. Therefore, some post process is inevitable. To make this step easier，different analysis results are attributed with different prefix of the line. For example, the mean d2min value in the log file starts with ' d|', so that we can get all the mean values using `grep` command as following
```bash
grep ' d|' out.txt > d2min.txt
```

2. A more complex situation is that some results include more than one line in the log file, e.g., the D2min distribution contains 400 lines after the title line, and the title line starts with ' d1|' which can be used as an indicator. Thus we need the following `grep` command to get all of those, and then split it into seperate files for each frame.
```bash
grep -A 400 ' d1|' out.txt > dis.txt
csplit -f my_ -n 3 dis.txt '/--/' '{*}' -s -k
```

3. There is an even more complex situation actually: the results include an uncertain number of lines. The starting and ending lines of the data have special starting strings, such as ' c1|', ' c2|'. Therefore, we can use `awk` to get all the related lines to a single file, and then use `csplit` command to devide it into independent files.
```bash
awk '/ c1\|/  {inblock = 1}  / c2\|/ {inblock = 0}  inblock == 1'  1mi_5e8_dc.log > dc.txt
csplit -f my_ -n 3 dc.txt '/ c1|/' '{*}' -s -k
```

