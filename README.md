# pbGraphs.ts
pbGraphs.ts is a graph library for TypeScript.

## Requirements
TypeScript 3.4+. The resulting JavaScript 5+ works in browsers, NodeJS and Java's script engine.

## Dependencies
There are no dependencies.

## Getting Started
It is easy to get started. Simply include the sources in your project. Check out test.ts in the src directory for an example.

## Avaialble in 12 different programming languages
pbGraphs-java has been implemented with [progsbase](https://www.progsbase.com), so the library is available for 12 different programming languages.

Try [all functions online](https://repo.progsbase.com/repoviewer/no.inductive.libraries/DirectedGraphs/0.1.14/) and see the implementations in different languages.

## Functions

### Graph Equality
```
function DirectedGraphsEqual(a : DirectedGraph, b : DirectedGraph) : boolean
```

Returns true of the two directed graphs are equal.

### Graph Components
```
function GetGraphComponents(g : DirectedGraph, componentMembership : NumberArrayReference) : boolean
```

Places the list of strongly connected components in componentMembership. Returns false if graphs does not quality as undirected.

### Topological Sort
```
function TopologicalSort(g : DirectedGraph, list : NumberArrayReference) : boolean
```

Places the topological sort of the graph in `list`.

### Searches

* Depth-first Search
```
function DepthFirstSearch(g : DirectedGraph, start : number, list : NumberArrayReference) : void
```

Places the depth-first sort ordering of the graph in `list`.

* Breadth-first Search
```
function BreadthFirstSearch(g : DirectedGraph, start : number, list : NumberArrayReference) : void
```

Places the breadth-first sort ordering of the graph in `list`.

### Shortest Paths

* Dijkstra's Algorithm
```
function DijkstrasAlgorithm(g : DirectedGraph, src : number, dist : NumberArrayReference, distSet : BooleanArrayReference, prev : NumberArrayReference) : void
```

Performs Dijkstra's algorithm on the graph `g` from `src`. Whether nodes are reachable is placed in `distSet`, the shortest distances in `dist`, and the previous node in the shortest paths in `prev`.

* Bellman-Ford Algorithm
```
function BellmanFordAlgorithm(g : DirectedGraph, src : number, dist : NumberArrayReference, distSet : BooleanArrayReference, prev : NumberArrayReference) : boolean
```

Performs the Bellman-Ford algorithm on the graph `g` from `src`. Whether nodes are reachable is placed in `distSet`, the shortest distances in `dist`, and the previous node in the shortest paths in `prev`.


* Floyd-Warshall Algorithm
```
function FloydWarshallAlgorithm(g : DirectedGraph, distances : Distances) : boolean
```

Performs the Floyd-Warshall algorithm on the graph `g`. The shortest distances between each pair of nodes are placed in `distances`.

### Minimum Spanning Trees

* Prim's Algorithm
```
function PrimsAlgorithm(g : DirectedGraph, forest : Forest) : boolean
```

Performs the Prim's algorithm on the graph `g`. All minimum spanning trees of the graph are placed in `forest`. Returns false if graphs does not quality as undirected.


* Kruskal's Algorithm
```
function KruskalsAlgorithm(g : DirectedGraph, forest : Forest) : boolean
```

Performs the Kruskal's algorithm on the graph `g`. All minimum spanning trees of the graph are placed in `forest`. Returns false if graphs does not quality as undirected.

### Cycles

* Cycle Detection
```
function DirectedGraphContainsCycleDFS(g : DirectedGraph) : boolean
```

Return true if there are cycles in the graph `g`.

* Cycle Counting
```
function DirectedGraphCountCyclesDFS(g : DirectedGraph) : number
```

Return the number of cycles in the graph `g`.

* Get All Cyles
```
function DirectedGraphGetCyclesDFS(g : DirectedGraph) : Cycle []
```

Returns the list of cycles in the graph `g`.