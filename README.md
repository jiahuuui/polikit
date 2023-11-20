# Polikit
## A polyhedral analysis toolkit

This package is originally developed for polyhedral analysis of amorphous structures. Now it has a module for topological constrain analysis also.

## Usage

'a.exe -f abc.xyz -p 1 -r 2.35 -o pt'
'-p' can be 1 or 0, decides whether periodic boundary condition will be applied.
'-r' is a float number, indicate the cutoff distance.
'-o' gives computing options, 'p' means polyhedral analysis, 't' means tct analysis.

'a.exe -d dumpfiles/*.xyz 20 -p 1 -r 2.35 -o pt'
'-d' gives a directory name that contains .xyz files, '*.xyz' is here to specify the file type. 20 is the frame invertal for dynamic comparison.

### Sep. 2023 -

2. Perform Delaunay triangulation on the system so that cavity analysis can be performed. Other methods: Gaussian density analysis, brutal force cutoff analysis.
  - Delauney triangulation can be achieved in two distinct ways: 1. from an initial small unit tetrahedron, insert point from outside of the convex hull once a time untill all points are inserted. 2. from a big enough initial tetrahedron that include all points inside it, insert point from inside of the tetrahedron once a time untill all points are inserted.
  - At PBC, atoms at outlayer need to be dealed seperately.
3. Perform ring statistics analysis. However the simplist breadth-first algorithm is super computing power expensive.

## Steps to perform a ring analysis

- find shortest path around a center atom
- check if an atom is in a pair of prime-mid-node
-
### 2/11/23
1. Now static analysis can be performed, but dynamic analysis still can't be performed. a buffer zone module is needed to store parameter needed to compare between frames(could use binary files).
~~2. How to read all files in a directory?~~

### 11/11/23
~~1. In dynamic mood, the fname should be changed each step, so that when calling the poly/tct analysis, the correct file would be opened. Then before comparison is actually performed, the numbers can be stored to a linked list. Other options are: save the data to be compared as binary file; or construct big enough array according to the interval read from keyboard.~~
2. Performing RDF analysis need exceptionally large cutoff, so the capacity of the bin should be increased. Though we have a bin capacity can be modified at this moment, error is that deallocate such a bincell is not allowed. We need to figure out other solutions.
3. When using linked list to store the computed ln/tct lists, the list size must be allocated after initialization according to the atom number. How should this be realized?
4. TCT analysis has not been tested yet.
5. read .dump format files.

<!-- 
### Heading
# H1
## H2
### H3
### Bold
**bold text**
### Italic
*italicized text*
### Blockquote
> blockquote
### Ordered List
1. First item
2. Second item
3. Third item
### Unordered List
- First item
- Second item
- Third item
### Code
`code`
### Horizontal Rule
---
### Link
[Markdown Guide](https://www.markdownguide.org)
### Image
![alt text](https://www.markdownguide.org/assets/images/tux.png)
### Table
| Syntax | Description |
| ----------- | ----------- |
| Header | Title |
| Paragraph | Text |
### Fenced Code Block
```fortran
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
```
### Strikethrough
~~The world is flat.~~
-->
