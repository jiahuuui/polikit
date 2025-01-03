
## 1. Topological Constraint Theory

Performing TCT analysis requires 2 extra parameters. First parameter is the frame range number during which the evaluation of standard deviation of bond length and bond angle is performed. Second is the threshold value beyond which we consider the constraint is inactive. TCT analysis requires intensive dumping frequency (short time seperation between frames) to get more reliable data.

**TCT analysis is performed in 3 steps:**

### 1.1. Construct container

If not initialized, two containers are created to store bond lengths and bond angles respectively. At each frame, the data from the earliest frame is removed and data from the newest frame is pushed into the container. If the TCT analysis is performed, then the memory won't be freed from frame to frame, but the data would be updated.

### 1.2. Clear space

Before storing data into the containters, the continuity of data should be checked first. This is because for each center atom, we prepared space for 10 sets of bond lengths and 30 sets of bond angles. If a bond breaks at one moment, the corresponding space should be freed as soon as possible, otherwise the space can run out faster and there would be no space for new data.

Unless the analysis has just started (current frame < frame range number), a data set should be removed as long as there is no new data at current frame. Because that means a broken bond. But when the analysis has just begun, the current frame number should be considered in the analysis.

### 1.3. Push in data

The data container has a header list and its corresponding data list. For example, for one center atom, the header list has a capacity of 10, which means it can have 10 sets of data in maximum. The header is pushed in order of ID.

When pushing new data into the container, itshould be pushed into the correct set. There can either be a found header, or a header to be inserted. It can be inserted to half the header list, or at the end of the list. After inserting the header, the new data would be pushed to the data set.

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
