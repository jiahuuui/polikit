## 2. Ring statistics analysis

When doing ring statistic analysis, the most important thing is the definition of the ring. This is because it determines how the analysis should be performed and affects the results significantly. There are some definitions that have been used:

**King’s criteria** - given two adjacent nodes (atoms), a ring is the shortest path from one node (atom) to the other.

**Shortest path criteria** - given a center node (atom) and two of its neighbors, a ring is the shortest path from one neighbor to the other without passing the center node (atom).

**Primitive ring criteria**^[https://doi.org/10.1016/0022-3093(90)90686-G]^[https://doi.org/10.1016/S0927-0256(01)00256-7] - a primitive ring is that the two paths between any node (atom) on the ring and its corresponding prime mid-node (the furthest node on the ring) are both shortest paths. In practice, larger rings can be determined.

**Strong ring criteria**^[https://doi.org/10.1016/0022-3093(91)90145-V] - a ring that can not be decomposed as a sum of any number of smaller rings is a strong ring. Definition of strong ring is  tricky, e.g., it is intuitional that Ring 1 in is a strong ring, and Ring 3 is not, but by flipping point A and B, Ring 2 and Ring 3 can be swapped, even Ring 1 can be regarded as topologically encloses Ring 2 and Ring 3. Hence, it is hard to understand how to identify a strong ring.

In this tool, the Primitive ring criteria is applied. The analysis needs a parameter to limit the maximum size of the ring, and also limit the memory cost. It is performed in following steps.

### 2.1 Simple ring statistics check

This method uses a simple and straightforward method to identify the primitive rings. The shortest paths list is created, the rings are formed and then checked if any shortcut can be found. If not, they are considered primitive. There are still two problems: the code gets slower because the ring list is getting larger; many repeat analysis is still performed.

#### 2.1.1 Create shortest path array

First, given a center atom, and the already constructed neighbor list, a breadth-first algorithm is used to create a shortest path array of the center atom. It stores all the possible shortest paths begin from the center to the allowed furthest atoms.

#### 2.1.2 Construct rings

In this step two types of rings are checked, i.e., even ring and odd ring. For an even ring, if an atom, namely *A*, appears twice in the same column (*n*) of the path list array, it is an even ring. For an odd ring, if atom *A* appears in column *n* and *n-1*, and correspondingly an atom *B* appears in the same row but inverse columns, it is an odd ring.

According to the primitive ring criteria, it is necessary that between a pair of nodes there must be two shortest paths. Therefore, from the shortest path array, the atoms that appear for more than once are considered end nodes of rings. All the rings are stored in a list for checking in next step.

#### 2.1.3 Remove non-primitive rings

For each ring found in last step, a shortest path search is performed on each pair of farthest nodes. This is to check if there is shortcut other than the two ring branches between them. If any shortcut is found, the ring is removed from the ring list, so that the rings left in the list are all primitive rings at the end of analysis.

The primitive ring search will be performed over all the atoms to find all the primitive rings.

### 2.2
<!--
There are three steps to find all the rings around a given node:
1. Build the shortest path array (SPA) of the center node.
2. Build the visibility array (VA) and insert the information of rings already found.
3. Search for new rings in the shortest path array and modify the visibility array-->

#### 2.2.1 Create shortest path array

The shortest path array (SPA) contains all the possible shortest paths starts from the
center node. Basically, it is an array that contains all the breadth-first search
results. It guarantees that all the rings identified are formed with two
shortest paths.

#### 2.2.2 Modify visibility array

Visibility array (VA) determines whether two shortest paths are ‘visible’ to each other.
At the beginning, most shortest paths are visible to each other, except the path and itself. Thus the VA is a matrix of all `.true.` but the diagonal elements are `.false.`. Given a ring, the VA can be modified in such a way that, the two branches of the ring and all the subsequent branches are 'invisible' to each other, which means all the corresponding elements are set to `.false.`.

VA is useful in two ways:

1. It can avoid repeat identifying of rings around a center atom. Before the analysis about the center atom, this is achieved by modifying the VA. The rings that contains the center atom are found, and the VA is modified according to the rings already found.

2. It can avoid further identifying of rings when a ring is already identified. During the analysis, the branches and subsequent branches are set 'invisible' to each other so that further check is avoided.

In practice, the way the VA is modified is also affected by what kind of rings we are looking for.

#### 2.2.3 Construct rings

If two shortest paths meet at the same node, a ring is identified. Then we
can set the two paths (or their upstream paths) invisible to each other to limit the
identifying of new rings. Ideally, at some point the visibility array will be fulfilled
with `.false.`, which means none of the shortest path is visible to any of them, and
we can stop searching for this center node.

The SPA can be large (~10^4 or ~10^5 shortest paths), I tried to dynamically modify
the SPA and VA so that the size of the arrays can be limited. However, in such a
way step 2 and 3 are merged and it is impossible to initialize the VA with known
rings, i.e., can not avoid repetition.

#### 2.2.4 Remove non-primitive rings
