# Polikit
## A polyhedral analysis toolkit

This package is originally developed for polyhedral analysis of amorphous structures. Now it has a module for analysis related to topological constraint theory also. It takes the dumped atomic configurations as input, and is able to perform static or dynamic analysis on the configurations.

## Usage
For static analysis:

**`a.exe -f abc.xyz -p 1 -r 2.35 -o pt`**

'-p' can be 1 or 0, decides whether periodic boundary condition will be applied.
'-r' is a float number, indicate the cutoff distance.
'-o' gives computing options, 'p' means polyhedral analysis, 't' means tct analysis.

For dynamic analysis:

**`a.exe -d dumpfiles 20 -p 1 -r 2.35 -o pt`**

'-d' gives a directory name that contains .xyz files. 20 is the frame invertal for dynamic comparison. Other parameters work in the same way as in static analysis.

### About TCT
Performing TCT analysis requires 2 extra parameters, and requires intensive dumping frequency (short frame interval) to get meaningful data. First parameter is the frame number range that the evaluation of standard deviation of bond length and bond angle is performed. Second is the threshold value beyond which we consider the constraint is inactive. These two values at this moment can only be configured in the source code (at `src/tct.f90`).

### Plans
1. Performing RDF analysis need exceptionally large cutoff, so the capacity of the bin should be increased. Though we have a bin capacity can be modified at this moment, error is that deallocate such a bincell is not allowed. We need to figure out other solutions.
3. read .dump format files.
4. find good way to dump atomic properties.
5. decide how the comparison data is exported.
6. l.133 main.f90

#### Delauney triangulation

Perform Delaunay triangulation on the system so that cavity analysis can be performed. Other methods: Gaussian density analysis, brutal force cutoff analysis.
  - Delauney triangulation can be achieved in two distinct ways: 1. from an initial small unit tetrahedron, insert point from outside of the convex hull once a time untill all points are inserted. 2. from a big enough initial tetrahedron that include all points inside it, insert point from inside of the tetrahedron once a time untill all points are inserted.
  - At PBC, atoms at outlayer need to be dealed seperately.
#### Ring statistics analysis

However the simplist breadth-first algorithm is super computing power expensive.

- find shortest path around a center atom
- check if an atom is in a pair of prime-mid-node

#### to Kai
Hello Kai,

Hope this email finds your well!

It has been nearly a year since I physically left our lab. I worked intensively with Antti this year to summarize my work and get ready for the doctoral defence. Hopefully I can get prepared for it in near future.

Without your introduction I would never have the chance to join our lab, I really feel grateful for it. In past years your attitude towards research always inspires me. Regarding my own research, we finished the project with many influential publications, including one publication on Advanced Materials and one on Acta Materialia which I made major contribution.

I found during my phd that my interest is still scientific research. More specifically, interested in topics include radiation effect in semiconductor materials such as gallium oxide. Therefore, I am wondering if it is possible to ask for a recommendation letter from you.

Happy holidays,
Jiahui


<!-- 
### Heading
# H1
## H2
### H3
#### Bold
**bold text**
##### Italic
*italicized text*
###### Blockquote
> blockquote
### Ordered List
1. First item
2. Second item
3. Third item
### Unordered List
- First item
- Second item
### highlight
=this is some highlight=
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
```java
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
```
### Strikethrough
~~The world is flat.~~
 -->
