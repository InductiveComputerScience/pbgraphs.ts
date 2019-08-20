
export class Distance{
	to : Target [];
}
export class Distances{
	fromx : Distance [];
}
export class Target{
	length : number;
	lengthSet : boolean;
	next : number;
}
export class Cycle{
	nodeNrs : number [];
}
export class DirectedGraph{
	nodes : Nodex [];
}
export class DirectedGraphMatrix{
	c : DirectedGraphMatrixColumn [];
}
export class DirectedGraphMatrixColumn{
	r : number [];
}
export class DirectedGraphReference{
	graph : DirectedGraph;
}
export class Edge{
	nodeNr : number;
	weight : number;
}
export class Nodex{
	edge : Edge [];
}
export class BooleanArrayReference{
	booleanArray : boolean [];
}
export class BooleanReference{
	booleanValue : boolean;
}
export class CharacterReference{
	characterValue : string;
}
export class NumberArrayReference{
	numberArray : number [];
}
export class NumberReference{
	numberValue : number;
}
export class StringArrayReference{
	stringArray : StringReference [];
}
export class StringReference{
	stringx : string [];
}
export class PriorityQueueBTNumbers{
	heap : DynamicArrayNumbers;
}
export class PriorityQueueBTNumKeyValue{
	heapKey : DynamicArrayNumbers;
	heapValue : DynamicArrayNumbers;
}
export class LinkedListNodeNumbers{
	next : LinkedListNodeNumbers;
	end : boolean;
	value : number;
}
export class LinkedListNumbers{
	first : LinkedListNodeNumbers;
	last : LinkedListNodeNumbers;
}
export class DynamicArrayNumbers{
	array : number [];
	length : number;
}
export class Forest{
	trees : Tree [];
}
export class Tree{
	branches : Tree [];
	label : number;
}
	export function DepthFirstSearch(g : DirectedGraph, start : number, list : NumberArrayReference) : void{
		var visited : boolean [];
		var ll : LinkedListNumbers;

		visited = CreateBooleanArray(g.nodes.length, false);
		ll = CreateLinkedListNumbers();

		DepthFirstSearchRecursive(g, g.nodes[start], start, visited, ll);

		list.numberArray = LinkedListNumbersToArray(ll);
		FreeLinkedListNumbers(ll);
	}


	export function DepthFirstSearchRecursive(g : DirectedGraph, node : Nodex, nodeNr : number, visited : boolean [], list : LinkedListNumbers) : void{
		var i : number;
		var e : Edge;

		visited[nodeNr] = true;

		LinkedListAddNumber(list, nodeNr);

		for(i = 0; i < node.edge.length; i = i + 1){
			e = node.edge[i];
			if(!visited[e.nodeNr]){
				DepthFirstSearchRecursive(g, g.nodes[e.nodeNr], e.nodeNr, visited, list);
			}
		}
	}


	export function BreadthFirstSearch(g : DirectedGraph, start : number, list : NumberArrayReference) : void{
		var visited : boolean [];
		var i : number, front : number, v : number, length : number;
		var e : Edge;
		var n : Nodex;
		var da : DynamicArrayNumbers;

		da = CreateDynamicArrayNumbers();
		visited = CreateBooleanArray(g.nodes.length, false);
		length = 0;
		front = 0;

		visited[start] = true;

		DynamicArrayAddNumber(da, start);
		length = length + 1;

		for(; front != length; ){
			v = DynamicArrayNumbersIndex(da, front);
			front = front + 1;

			n = g.nodes[v];

			for(i = 0; i < n.edge.length; i = i + 1){
				e = n.edge[i];
				if(!visited[e.nodeNr]){
					visited[e.nodeNr] = true;

					DynamicArrayAddNumber(da, e.nodeNr);
					length = length + 1;
				}
			}
		}

		list.numberArray = DynamicArrayNumbersToArray(da);
		FreeDynamicArrayNumbers(da);
	}


	export function PrimsAlgorithmNoQueue(g : DirectedGraph, forest : Forest) : boolean{
		var valid : boolean, found : boolean, minimumSet : boolean;
		var inMST : boolean [];
		var i : number, j : number, root : number, minimum : number, minimumTarget : number, minimumSource : number, nodesCompleted : number;
		var node : Nodex;
		var edge : Edge;
		var linkedListTrees : LinkedListNumbers [];
		var roots : LinkedListNumbers;

		valid = DirectedGraphIsValid(g) && IsUndirected(g);

		if(valid){
			inMST = CreateBooleanArray(g.nodes.length, false);
			nodesCompleted = 0;
			linkedListTrees = CreateLinkedListNumbersArray(g.nodes.length);
			roots = CreateLinkedListNumbers();

			for(; nodesCompleted < g.nodes.length; ){

				/* Find a node not in an MST*/
				found = false;
				root = 0;
				for(i = 0; i < g.nodes.length && !found; i = i + 1){
					if(!inMST[i]){
						root = i;
						found = true;
					}
				}

				LinkedListAddNumber(roots, root);

				inMST[root] = true;
				nodesCompleted = nodesCompleted + 1;

				found = true;
				for(; found; ){
					/* Find minimum edge going out from existing tree.*/
					minimum = 0;
					minimumSet = false;
					minimumTarget = 0;
					minimumSource = 0;
					for(i = 0; i < g.nodes.length; i = i + 1){
						if(inMST[i]){
							node = g.nodes[i];
							for(j = 0; j < node.edge.length; j = j + 1){
								edge = node.edge[j];
								if(!inMST[edge.nodeNr]){
									if(!minimumSet){
										minimum = edge.weight;
										minimumTarget = edge.nodeNr;
										minimumSource = i;
										minimumSet = true;
									}else if(edge.weight < minimum){
										minimum = edge.weight;
										minimumTarget = edge.nodeNr;
										minimumSource = i;
									}
								}
							}
						}
					}

					/* Add edge to tree.*/
					if(minimumSet){
						LinkedListAddNumber(linkedListTrees[minimumSource], minimumTarget);
						inMST[minimumTarget] = true;
						nodesCompleted = nodesCompleted + 1;
						found = true;
					}else{
						found = false;
					}
				}
			}

			ConvertLinkedListTreesToForest(forest, roots, linkedListTrees);

			/* Free memory.*/
			inMST = undefined;
			FreeLinkedListNumbersArray(linkedListTrees);
		}

		return valid;
	}


	export function PrimsAlgorithm(g : DirectedGraph, forest : Forest) : boolean{
		var valid : boolean, found : boolean, minimumSet : boolean, empty : boolean;
		var inMST : boolean [], minimumEdgeSet : boolean [];
		var i : number, root : number, minimumTarget : number, minimumSource : number, nodesCompleted : number;
		var minimumEdges : number [], minimumSources : number [];
		var q : PriorityQueueBTNumKeyValue;
		var targetReference : NumberReference, weightReference : NumberReference;
		var linkedListTrees : LinkedListNumbers [];
		var roots : LinkedListNumbers;

		valid = DirectedGraphIsValid(g) && IsUndirected(g);

		if(valid){
			q = CreatePriorityQueueBTNumKeyValue();
			targetReference = new NumberReference();
			weightReference = new NumberReference();
			inMST = CreateBooleanArray(g.nodes.length, false);
			minimumEdgeSet = CreateBooleanArray(g.nodes.length, false);
			minimumEdges = CreateNumberArray(g.nodes.length, 0);
			minimumSources = CreateNumberArray(g.nodes.length, 0);
			linkedListTrees = CreateLinkedListNumbersArray(g.nodes.length);
			roots = CreateLinkedListNumbers();
			nodesCompleted = 0;

			for(; nodesCompleted < g.nodes.length; ){
				/* Find a node not in an MST*/
				found = false;
				root = 0;
				for(i = 0; i < g.nodes.length && !found; i = i + 1){
					if(!inMST[i]){
						root = i;
						found = true;
					}
				}

				/* Record tree root.*/
				LinkedListAddNumber(roots, root);
				inMST[root] = true;
				nodesCompleted = nodesCompleted + 1;

				/* Add all outgoing edges to priority queue*/
				AddOutgoingEdgesToPriorityQueue(g.nodes[root], minimumEdgeSet, root, minimumEdges, minimumSources, q);

				/* Expand tree one vertex at a time.*/
				found = true;
				for(; found; ){
					/* Find minimum edge going out from existing tree using queue.*/
					minimumSet = false;
					empty = false;
					for(; !minimumSet && !empty; ){
						empty = !PopPriorityQueueBTNumKeyValue(q, weightReference, targetReference);
						if(!empty && !inMST[targetReference.numberValue]){
							minimumSet = true;
						}
					}

					if(minimumSet){
						/* Add edge to tree.*/
						minimumTarget = targetReference.numberValue;
						minimumSource = minimumSources[minimumTarget];

						LinkedListAddNumber(linkedListTrees[minimumSource], minimumTarget);
						inMST[minimumTarget] = true;
						nodesCompleted = nodesCompleted + 1;
						found = true;

						/* Add all outgoing edges to priority queue.*/
						AddOutgoingEdgesToPriorityQueue(g.nodes[minimumTarget], minimumEdgeSet, minimumTarget, minimumEdges, minimumSources, q);
					}else{
						found = false;
					}
				}
			}

			ConvertLinkedListTreesToForest(forest, roots, linkedListTrees);

			/* Free memory.*/
			FreePriorityQueueBTNumKeyValue(q);
			targetReference = undefined;
			weightReference = undefined;
			inMST = undefined;
			minimumEdgeSet = undefined;
			FreeLinkedListNumbersArray(linkedListTrees);
			minimumEdges = undefined;
			minimumSources = undefined;
		}

		return valid;
	}


	export function AddOutgoingEdgesToPriorityQueue(node : Nodex, minimumEdgeSet : boolean [], source : number, minimumEdges : number [], minimumSources : number [], q : PriorityQueueBTNumKeyValue) : void{
		var i : number, target : number;
		var edge : Edge;

		for(i = 0; i < node.edge.length; i = i + 1){
			edge = node.edge[i];
			target = edge.nodeNr;
			InsertIntoPriorityQueueBTNumKeyValue(q, 1/edge.weight, target);
			if(!minimumEdgeSet[target]){
				minimumEdges[target] = edge.weight;
				minimumSources[target] = source;
				minimumEdgeSet[target] = true;
			}else if(minimumEdges[target] > edge.weight){
				minimumEdges[target] = edge.weight;
				minimumSources[target] = source;
			}
		}
	}


	export function KruskalsAlgorithm(g : DirectedGraph, forest : Forest) : boolean{
		var valid : boolean;
		var q : PriorityQueueBTNumKeyValue;
		var node : Nodex;
		var edge : Edge;
		var i : number, j : number, edgeNr : number, source : number, target : number, replace : number, replaceWith : number, candidate : number, trees : number, treeNr : number;
		var sources : DynamicArrayNumbers, targets : DynamicArrayNumbers;
		var edges : DynamicArrayNumbers;
		var memberOfTree : number [];
		var roots : boolean [];
		var edgeNrReference : NumberReference, weightReference : NumberReference;
		var tree : Tree;

		valid = DirectedGraphIsValid(g) && IsUndirected(g);

		if(valid){
			sources = CreateDynamicArrayNumbers();
			targets = CreateDynamicArrayNumbers();
			edges = CreateDynamicArrayNumbers();
			edgeNrReference = new NumberReference();
			weightReference = new NumberReference();
			roots = CreateBooleanArray(g.nodes.length, false);
			memberOfTree = new Array<number>(g.nodes.length);
			for(i = 0; i < g.nodes.length; i = i + 1){
				memberOfTree[i] = i;
			}

			q = CreatePriorityQueueBTNumKeyValue();

			/* Add all edges to a priority queue.*/
			edgeNr = 0;
			for(i = 0; i < g.nodes.length; i = i + 1){
				node = g.nodes[i];
				for(j = 0; j < node.edge.length; j = j + 1){
					edge = node.edge[j];
					InsertIntoPriorityQueueBTNumKeyValue(q, 1/edge.weight, edgeNr);
					DynamicArrayAddNumber(sources, i);
					DynamicArrayAddNumber(targets, edge.nodeNr);

					edgeNr = edgeNr + 1;
				}
			}

			for(; !IsEmptyPriorityQueueBTNumKeyValue(q); ){
				PopPriorityQueueBTNumKeyValue(q, weightReference, edgeNrReference);

				source = DynamicArrayNumbersIndex(sources, edgeNrReference.numberValue);
				target = DynamicArrayNumbersIndex(targets, edgeNrReference.numberValue);

				if(memberOfTree[source] != memberOfTree[target]){
					replace = memberOfTree[target];
					replaceWith = memberOfTree[source];

					for(i = 0; i < g.nodes.length; i = i + 1){
						if(memberOfTree[i] == replace){
							memberOfTree[i] = replaceWith;
						}
					}

					DynamicArrayAddNumber(edges, edgeNrReference.numberValue);
				}
			}

			/* Built forest.*/
			trees = 0;
			for(i = 0; i < g.nodes.length; i = i + 1){
				candidate = memberOfTree[i];
				if(!roots[candidate]){
					trees = trees + 1;
					roots[candidate] = true;
				}
			}
			forest.trees = new Array<Tree>(trees);
			treeNr = 0;
			for(i = 0; i < g.nodes.length; i = i + 1){
				if(roots[i]){
					tree = CreateTreeFromEdgeCollection(i, i, edges, sources, targets);

					forest.trees[treeNr] = tree;
					treeNr = treeNr + 1;
				}
			}

			/* Free memory.*/
			FreePriorityQueueBTNumKeyValue(q);
			FreeDynamicArrayNumbers(sources);
			FreeDynamicArrayNumbers(targets);
			edgeNrReference = undefined;
			weightReference = undefined;
		}

		return valid;
	}


	export function CreateTreeFromEdgeCollection(root : number, parent : number, edges : DynamicArrayNumbers, sources : DynamicArrayNumbers, targets : DynamicArrayNumbers) : Tree{
		var tree : Tree;
		var i : number, edgeNr : number, source : number, target : number, size : number;
		var branches : LinkedListNumbers;
		var node : LinkedListNodeNumbers;

		tree = new Tree();
		tree.label = root;
		branches = CreateLinkedListNumbers();

		size = 0;
		for(i = 0; i < edges.length; i = i + 1){
			edgeNr = DynamicArrayNumbersIndex(edges, i);

			source = DynamicArrayNumbersIndex(sources, edgeNr);
			target = DynamicArrayNumbersIndex(targets, edgeNr);

			if(source == root && target != parent){
				LinkedListAddNumber(branches, target);
				size = size + 1;
			}else if(target == root && source != parent){
				LinkedListAddNumber(branches, source);
				size = size + 1;
			}
		}

		tree.branches = new Array<Tree>(size);
		node = branches.first;
		for(i = 0; i < size; i = i + 1){
			tree.branches[i] = CreateTreeFromEdgeCollection(node.value, root, edges, sources, targets);

			node = node.next;
		}

		return tree;
	}


	export function DijkstrasAlgorithm(g : DirectedGraph, src : number, dist : NumberArrayReference, distSet : BooleanArrayReference, prev : NumberArrayReference) : void{
		var nodeDone : boolean [];
		var i : number, j : number, v : number, nodes : number, u : number, edges : number, alt : number;
		var edge : Edge;

		nodes = g.nodes.length;

		distSet.booleanArray = CreateBooleanArray(nodes, false);
		nodeDone = CreateBooleanArray(nodes, false);
		dist.numberArray = CreateNumberArray(nodes, 0);
		prev.numberArray = CreateNumberArray(nodes, 0);

		dist.numberArray[src] = 0;
		distSet.booleanArray[src] = true;

		for(i = 0; i < nodes; i = i + 1){
			/* Get node with lowest distance*/
			u = ListFindLowestSetAndIncluded(dist.numberArray, distSet.booleanArray, nodeDone);

			/* Mark node as done*/
			nodeDone[u] = true;

			edges = GetEdgesForNodeFromDirectedGraph(g, u);
			for(j = 0; j < edges; j = j + 1){
				edge = GetEdgeFromDirectedGraph(g, u, j);

				if(!nodeDone[edge.nodeNr]){
					v = edge.nodeNr;
					alt = dist.numberArray[u] + edge.weight;
					if(!distSet.booleanArray[v]){
						dist.numberArray[v] = alt;
						distSet.booleanArray[v] = true;
						prev.numberArray[v] = u;
					}else if(alt < dist.numberArray[v]){
						dist.numberArray[v] = alt;
						prev.numberArray[v] = u;
					}
				}
			}
		}

		distSet = undefined;
		nodeDone = undefined;
	}


	export function ListFindLowestSetAndIncluded(list : number [], setx : boolean [], exclude : boolean []) : number{
		var i : number, nodes : number, lowest : number, u : number;
		var lowestSet : boolean;

		nodes = list.length;
		lowest = 0;
		u = 0;

		lowestSet = false;
		for(i = 0; i < nodes; i = i + 1){
			if(!exclude[i] && setx[i]){
				if(!lowestSet){
					lowest = list[i];
					u = i;
					lowestSet = true;
				}else if(list[i] < lowest){
					lowest = list[i];
					u = i;
				}
			}
		}

		return u;
	}


	export function DijkstrasAlgorithmDestinationOnly(g : DirectedGraph, src : number, dest : number, path : NumberArrayReference, distance : NumberReference) : boolean{
		var distances : NumberArrayReference, previous : NumberArrayReference;
		var distanceSet : BooleanArrayReference;
		var found : boolean;

		distances = new NumberArrayReference();
		previous = new NumberArrayReference();
		distanceSet = new BooleanArrayReference();

		DijkstrasAlgorithm(g, src, distances, distanceSet, previous);

		found = distanceSet.booleanArray[dest];

		if(found){
			distance.numberValue = distances.numberArray[dest];

			ExtractForwardPathFromReverseList(src, dest, previous, path);
		}

		distances = undefined;
		previous = undefined;
		distanceSet = undefined;

		return found;
	}


	export function ExtractForwardPathFromReverseList(src : number, dest : number, previous : NumberArrayReference, path : NumberArrayReference) : void{
		var next : number, length : number, i : number;

		next = dest;
		for(length = 1; next != src; length = length + 1){
			next = previous.numberArray[next];
		}

		path.numberArray = CreateNumberArray(length, 0);

		next = dest;
		for(i = 0; i < length; i = i + 1){
			path.numberArray[length - i - 1] = next;
			next = previous.numberArray[next];
		}
	}


	export function BellmanFordAlgorithm(g : DirectedGraph, src : number, dist : NumberArrayReference, distSet : BooleanArrayReference, prev : NumberArrayReference) : boolean{
		var nodeDone : boolean [];
		var i : number, j : number, v : number, nodes : number, u : number, edges : number, w : number;
		var edge : Edge;
		var success : boolean;

		nodes = g.nodes.length;

		distSet.booleanArray = CreateBooleanArray(nodes, false);
		nodeDone = CreateBooleanArray(nodes, false);
		dist.numberArray = CreateNumberArray(nodes, 0);
		prev.numberArray = CreateNumberArray(nodes, 0);

		dist.numberArray[src] = 0;
		distSet.booleanArray[src] = true;

		for(i = 0; i < nodes - 1; i = i + 1){
			for(u = 0; u < nodes; u = u + 1){
				edges = GetEdgesForNodeFromDirectedGraph(g, u);
				for(j = 0; j < edges; j = j + 1){
					edge = GetEdgeFromDirectedGraph(g, u, j);

					v = edge.nodeNr;
					w = edge.weight;

					if(distSet.booleanArray[u]){
						if(!distSet.booleanArray[v]){
							dist.numberArray[v] = dist.numberArray[u] + w;
							distSet.booleanArray[v] = true;
							prev.numberArray[v] = u;
						}else if(dist.numberArray[u] + w < dist.numberArray[v]){
							dist.numberArray[v] = dist.numberArray[u] + w;
							prev.numberArray[v] = u;
						}
					}
				}
			}
		}

		success = true;
		for(u = 0; u < nodes; u = u + 1){
			edges = GetEdgesForNodeFromDirectedGraph(g, u);
			for(j = 0; j < edges; j = j + 1){
				edge = GetEdgeFromDirectedGraph(g, u, j);

				v = edge.nodeNr;
				w = edge.weight;

				if(dist.numberArray[u] + w < dist.numberArray[v]){
					success = false;
				}
			}
		}

		distSet = undefined;
		nodeDone = undefined;

		return success;
	}


	export function BellmanFordAlgorithmDestinationOnly(g : DirectedGraph, src : number, dest : number, path : NumberArrayReference, distance : NumberReference) : boolean{
		var distances : NumberArrayReference, previous : NumberArrayReference;
		var distanceSet : BooleanArrayReference;
		var found : boolean;

		distances = new NumberArrayReference();
		previous = new NumberArrayReference();
		distanceSet = new BooleanArrayReference();

		found = BellmanFordAlgorithm(g, src, distances, distanceSet, previous);

		if(found){
			found = distanceSet.booleanArray[dest];

			if(found){
				distance.numberValue = distances.numberArray[dest];

				ExtractForwardPathFromReverseList(src, dest, previous, path);
			}
		}

		distances = undefined;
		previous = undefined;
		distanceSet = undefined;

		return found;
	}


	export function FloydWarshallAlgorithm(g : DirectedGraph, distances : Distances) : boolean{
		var u : number, v : number, k : number, i : number, j : number;
		var n : Nodex;
		var e : Edge;
		var t : Target, ij : Target, ik : Target, kj : Target;
		var success : boolean;

		success = true;

		for(u = 0; u < g.nodes.length; u = u + 1){
			n = g.nodes[u];

			for(j = 0; j < n.edge.length; j = j + 1){
				e = n.edge[j];
				v = e.nodeNr;

				t = distances.fromx[u].to[v];
				t.length = e.weight;
				t.lengthSet = true;
				t.next = v;
			}
		}

		for(v = 0; v < g.nodes.length; v = v + 1){
			t = distances.fromx[v].to[v];
			t.length = 0;
			t.lengthSet = true;
			t.next = v;
		}

		for(k = 0; k < g.nodes.length && success; k = k + 1){
			for(i = 0; i < g.nodes.length && success; i = i + 1){
				for(j = 0; j < g.nodes.length && success; j = j + 1){
					ij = distances.fromx[i].to[j];
					ik = distances.fromx[i].to[k];
					kj = distances.fromx[k].to[j];

					if(!ij.lengthSet && ik.lengthSet && kj.lengthSet){
						ij.length = ik.length + kj.length;
						ij.lengthSet = true;
						ij.next = ik.next;
					}else if(ij.lengthSet && ik.lengthSet && kj.lengthSet){
						if(ij.length > ik.length + kj.length){
							ij.length = ik.length + kj.length;
							ij.next = ik.next;
						}
					}

					if(i == j){
						if(ij.length < 0){
							success = false;
						}
					}
				}
			}
		}

		return success;
	}


	export function CreateDistancesFloydWarshallAlgorithm(nodes : number) : Distances{
		var distances : Distances;
		var i : number, j : number;

		distances = new Distances();
		distances.fromx = new Array<Distance>(nodes);
		for(i = 0; i < distances.fromx.length; i = i + 1){
			distances.fromx[i] = new Distance();
			distances.fromx[i].to = new Array<Target>(distances.fromx.length);
			for(j = 0; j < distances.fromx.length; j = j + 1){
				distances.fromx[i].to[j] = new Target();
				distances.fromx[i].to[j].length = 0;
				distances.fromx[i].to[j].lengthSet = false;
				distances.fromx[i].to[j].next = 0;
			}
		}

		return distances;
	}


	export function GetPathFromDistances(distances : Distances, u : number, v : number) : number []{
		var path : number [];
		var length : number, next : number, i : number;
		var t : Target;

		t = distances.fromx[u].to[v];
		if(t.lengthSet){
			/* count*/
			length = 1;
			next = u;
			for(; next != v; ){
				next = distances.fromx[next].to[v].next;
				length = length + 1;
			}

			path = new Array<number>(length);

			/* set*/
			next = u;
			for(i = 0; i < length; i = i + 1){
				path[i] = next;
				next = distances.fromx[next].to[v].next;
			}
		}else{
			path = new Array<number>(0);
		}

		return path;
	}


	export function CreateEdge(nodeNr : number, weight : number) : Edge{
		var e : Edge;

		e = new Edge();

		e.nodeNr = nodeNr;
		e.weight = weight;

		return e;
	}


	export function DirectedGraphIsValid(g : DirectedGraph) : boolean{
		var valid : boolean;
		var i : number, j : number;
		var node : Nodex;
		var edge : Edge;

		valid = true;

		for(i = 0; i < g.nodes.length; i = i + 1){
			node = g.nodes[i];
			for(j = 0; j < node.edge.length; j = j + 1){
				edge = node.edge[j];
				if(IsInteger(edge.nodeNr)){
					if(edge.nodeNr >= 0 && edge.nodeNr < g.nodes.length){
					}else{
						valid = false;
					}
				}else{
					valid = false;
				}
			}
		}

		return valid;
	}


	export function DirectedGraphContainsCycleDFS(g : DirectedGraph) : boolean{
		var incoming : number [];
		var i : number, zeroIncomming : number;
		var cycle : boolean;

		cycle = false;
		incoming = DirectedGraphCountIncomingEdges(g);

		zeroIncomming = 0;
		for(i = 0; i < g.nodes.length && !cycle; i = i + 1){
			if(incoming[i] == 0){
				zeroIncomming = zeroIncomming + 1;

				cycle = DirectedGraphContainsCycleFromNodeDFS(g, i);
			}
		}

		incoming = undefined;

		if(g.nodes.length > 0 && zeroIncomming == 0){
			cycle = true;
		}

		return cycle;
	}


	export function DirectedGraphCountIncomingEdges(g : DirectedGraph) : number []{
		var incoming : number [];
		var i : number, j : number;
		var node : Nodex;
		var e : Edge;

		incoming = CreateNumberArray(g.nodes.length, 0);

		for(i = 0; i < g.nodes.length; i = i + 1){
			node = g.nodes[i];
			for(j = 0; j < node.edge.length; j = j + 1){
				e = node.edge[j];
				incoming[e.nodeNr] = incoming[e.nodeNr] + 1;
			}
		}

		return incoming;
	}


	export function DirectedGraphContainsCycleFromNodeDFS(g : DirectedGraph, nodeNr : number) : boolean{
		var visited : boolean [];
		var cycle : boolean;

		visited = CreateBooleanArray(g.nodes.length, false);

		cycle = DirectedGraphContainsCycleFromNodeDFSRecursive(g, nodeNr, visited);

		visited = undefined;

		return cycle;
	}


	export function DirectedGraphContainsCycleFromNodeDFSRecursive(g : DirectedGraph, nodeNr : number, visited : boolean []) : boolean{
		var i : number;
		var e : Edge;
		var cycle : boolean;
		var node : Nodex;

		cycle = false;
		node = g.nodes[nodeNr];

		for(i = 0; i < node.edge.length && !cycle; i = i + 1){
			e = node.edge[i];
			if(visited[e.nodeNr]){
				cycle = true;
			}else{
				visited[e.nodeNr] = true;
				cycle = DirectedGraphContainsCycleFromNodeDFSRecursive(g, e.nodeNr, visited);
				visited[e.nodeNr] = false;
			}
		}

		return cycle;
	}


	export function DirectedGraphCountCyclesDFS(g : DirectedGraph) : number{
		var incoming : number [];
		var i : number, zeroIncoming : number, cycleCount : number;

		cycleCount = 0;
		incoming = DirectedGraphCountIncomingEdges(g);

		zeroIncoming = 0;
		for(i = 0; i < g.nodes.length; i = i + 1){
			if(incoming[i] == 0){
				zeroIncoming = zeroIncoming + 1;

				cycleCount = cycleCount + DirectedGraphCountCyclesFromNodeDFS(g, i);
			}
		}

		if(g.nodes.length > 0 && zeroIncoming == 0){
			cycleCount = cycleCount + DirectedGraphCountCyclesFromNodeDFS(g, 0);
		}

		incoming = undefined;

		return cycleCount;
	}


	export function DirectedGraphCountCyclesFromNodeDFS(g : DirectedGraph, nodeNr : number) : number{
		var color : number [];
		var cycleCount : number;

		color = CreateNumberArray(g.nodes.length, 0);

		cycleCount = DirectedGraphCountCyclesFromNodeDFSRecursive(g, nodeNr, color);

		color = undefined;

		return cycleCount;
	}


	export function DirectedGraphCountCyclesFromNodeDFSRecursive(g : DirectedGraph, nodeNr : number, color : number []) : number{
		var i : number, cycleCount : number;
		var e : Edge;
		var node : Nodex;

		cycleCount = 0;
		node = g.nodes[nodeNr];

		color[nodeNr] = 1;

		for(i = 0; i < node.edge.length; i = i + 1){
			e = node.edge[i];

			if(color[e.nodeNr] != 2){
				if(color[e.nodeNr] == 1){
					cycleCount = cycleCount + 1;
				}else{
					cycleCount = cycleCount + DirectedGraphCountCyclesFromNodeDFSRecursive(g, e.nodeNr, color);
				}
			}
		}

		color[nodeNr] = 2;

		return cycleCount;
	}


	export function DirectedGraphGetCyclesDFS(g : DirectedGraph) : Cycle []{
		var cycleCount : number;
		var cycles : Cycle [];
		var incoming : number [];
		var i : number, zeroIncoming : number;
		var cycleNumber : NumberReference;

		cycleNumber = CreateNumberReference(0);
		cycleCount = DirectedGraphCountCyclesDFS(g);

		cycles = new Array<Cycle>(cycleCount);

		incoming = DirectedGraphCountIncomingEdges(g);

		zeroIncoming = 0;
		for(i = 0; i < g.nodes.length; i = i + 1){
			if(incoming[i] == 0){
				zeroIncoming = zeroIncoming + 1;

				DirectedGraphGetCyclesFromNodeDFS(g, i, cycles, cycleNumber);
			}
		}

		if(g.nodes.length > 0 && zeroIncoming == 0){
			DirectedGraphGetCyclesFromNodeDFS(g, 0, cycles, cycleNumber);
		}

		incoming = undefined;

		return cycles;
	}


	export function DirectedGraphGetCyclesFromNodeDFS(g : DirectedGraph, nodeNr : number, cycles : Cycle [], cycleNumber : NumberReference) : void{
		var color : number [], cycleMark : number [];
		var previous : number [];
		var previousLength : number;

		color = CreateNumberArray(g.nodes.length, 0);
		cycleMark = CreateNumberArray(g.nodes.length, 0);
		previous = CreateNumberArray(g.nodes.length, 0);
		previousLength = 0;

		DirectedGraphGetCyclesFromNodeDFSRecursive(g, nodeNr, cycleNumber, color, cycles, previous, previousLength);

		color = undefined;
		cycleMark = undefined;
	}


	export function DirectedGraphGetCyclesFromNodeDFSRecursive(g : DirectedGraph, nodeNr : number, cycleNumber : NumberReference, color : number [], cycles : Cycle [], previous : number [], previousLength : number) : void{
		var i : number, j : number, current : number, cycleLength : number, next : number;
		var e : Edge;
		var node : Nodex;
		var done : boolean;

		node = g.nodes[nodeNr];

		color[nodeNr] = 1;

		previous[previousLength] = nodeNr;

		for(i = 0; i < node.edge.length; i = i + 1){
			e = node.edge[i];
			if(color[e.nodeNr] != 2){
				if(color[e.nodeNr] == 1){
					/* Get cycle length*/
					cycleLength = 0;
					done = false;
					current = previousLength;
					for(; !done; ){
						cycleLength = cycleLength + 1;
						if(previous[current] == e.nodeNr){
							done = true;
						}
						current = current - 1;
					}

					/* Get cycle in order*/
					cycles[cycleNumber.numberValue] = new Cycle();
					cycles[cycleNumber.numberValue].nodeNrs = new Array<number>(cycleLength);
					for(j = 0; j < cycleLength; j = j + 1){
						next = previousLength - cycleLength + 1 + j;
						cycles[cycleNumber.numberValue].nodeNrs[j] = previous[next];
					}

					cycleNumber.numberValue = cycleNumber.numberValue + 1;
				}else{
					DirectedGraphGetCyclesFromNodeDFSRecursive(g, e.nodeNr, cycleNumber, color, cycles, previous, previousLength + 1);
				}
			}
		}

		color[nodeNr] = 2;
	}


	export function CreateDirectedGraph(nodes : number) : DirectedGraph{
		var directedGraph : DirectedGraph;
		var i : number;

		directedGraph = new DirectedGraph();
		directedGraph.nodes = new Array<Nodex>(nodes);

		for(i = 0; i < nodes; i = i + 1){
			directedGraph.nodes[i] = new Nodex();
		}

		return directedGraph;
	}


	export function CreateDirectedGraphFromMatrixForm(m : DirectedGraphMatrix) : DirectedGraph{
		var g : DirectedGraph;
		var nodes : number, i : number, j : number, order : number, edgeValue : number, edgeNr : number;

		nodes = GetNodesFromDirectedGraphMatrix(m);

		g = CreateDirectedGraph(nodes);

		for(i = 0; i < nodes; i = i + 1){
			order = GetNodeOrderFromMatrixForm(m, i);
			g.nodes[i].edge = new Array<Edge>(order);
			edgeNr = 0;
			for(j = 0; j < nodes; j = j + 1){
				edgeValue = GetDirectedGraphMatrixEdge(m, i, j);
				if(edgeValue != 0){
					g.nodes[i].edge[edgeNr] = CreateEdge(j, edgeValue);
					edgeNr = edgeNr + 1;
				}
			}
		}

		return g;
	}


	export function GetDirectedGraphMatrixEdge(m : DirectedGraphMatrix, nodeNr : number, edgeNr : number) : number{
		return m.c[nodeNr].r[edgeNr];
	}


	export function GetNodeOrderFromMatrixForm(m : DirectedGraphMatrix, nodeNr : number) : number{
		var nodes : number, i : number, order : number;

		nodes = GetNodesFromDirectedGraphMatrix(m);

		order = 0;
		for(i = 0; i < nodes; i = i + 1){
			order = order + m.c[nodeNr].r[i];
		}

		return order;
	}


	export function GetNodesFromDirectedGraphMatrix(m : DirectedGraphMatrix) : number{
		return m.c.length;
	}


	export function CreateDirectedGraphMatrixFromListForm(g : DirectedGraph) : DirectedGraphMatrix{
		var m : DirectedGraphMatrix;
		var i : number, j : number, nodes : number;
		var node : Nodex;
		var edge : Edge;

		nodes = g.nodes.length;
		m = CreateDirectedGraphMatrix(nodes);

		for(i = 0; i < nodes; i = i + 1){
			node = g.nodes[i];
			for(j = 0; j < node.edge.length; j = j + 1){
				edge = node.edge[j];
				m.c[i].r[edge.nodeNr] = edge.weight;
			}
		}

		return m;
	}


	export function CreateDirectedGraphMatrix(nodes : number) : DirectedGraphMatrix{
		var m : DirectedGraphMatrix;
		var i : number;
		m = new DirectedGraphMatrix();

		m.c = new Array<DirectedGraphMatrixColumn>(nodes);
		for(i = 0; i < nodes; i = i + 1){
			m.c[i] = new DirectedGraphMatrixColumn();
			m.c[i].r = CreateNumberArray(nodes, 0);
		}

		return m;
	}


	export function DirectedGraphsEqual(a : DirectedGraph, b : DirectedGraph) : boolean{
		var equal : boolean, found : boolean, done : boolean;
		var nodes : number, foundCount : number, i : number, j : number, k : number, edges : number;
		var edgeA : Edge, edgeB : Edge;

		equal = true;

		if(a.nodes.length == b.nodes.length){
			nodes = a.nodes.length;

			done = false;
			for(i = 0; i < nodes && !done; i = i + 1){
				if(GetEdgesForNodeFromDirectedGraph(a, i) == GetEdgesForNodeFromDirectedGraph(b, i)){
					edges = GetEdgesForNodeFromDirectedGraph(a, i);
					foundCount = 0;
					for(j = 0; j < edges && !done; j = j + 1){
						found = false;
						for(k = 0; k < edges && !found; k = k + 1){
							edgeA = GetEdgeFromDirectedGraph(a, i, j);
							edgeB = GetEdgeFromDirectedGraph(b, i, k);
							if(edgeA.nodeNr == edgeB.nodeNr){
								if(edgeA.weight == edgeB.weight){
									found = true;
								}
							}
						}
						if(found){
							foundCount = foundCount + 1;
						}else{
							equal = false;
							done = true;
						}
					}

					if(foundCount == edges){
					}else{
						equal = false;
						done = true;
					}
				}else{
					equal = false;
					done = true;
				}
			}
		}else{
			equal = false;
		}

		return equal;
	}


	export function GetEdgesForNodeFromDirectedGraph(g : DirectedGraph, nodeNr : number) : number{
		return g.nodes[nodeNr].edge.length;
	}


	export function GetEdgeFromDirectedGraph(g : DirectedGraph, nodeNr : number, edgeNr : number) : Edge{
		return g.nodes[nodeNr].edge[edgeNr];
	}


	export function DirectedGraphMatricesEqual(a : DirectedGraphMatrix, b : DirectedGraphMatrix) : boolean{
		var equal : boolean;
		var nodes : number, i : number, j : number;

		equal = true;

		if(GetNodesFromDirectedGraphMatrix(a) == GetNodesFromDirectedGraphMatrix(b)){
			nodes = GetNodesFromDirectedGraphMatrix(a);
			for(i = 0; i < nodes && equal; i = i + 1){
				for(j = 0; j < nodes && equal; j = j + 1){
					if(GetDirectedGraphMatrixEdge(a, i, j) == GetDirectedGraphMatrixEdge(b, i, j)){
					}else{
						equal = false;
					}
				}
			}
		}else{
			equal = false;
		}

		return equal;
	}


	export function IsUndirected(g : DirectedGraph) : boolean{
		var undirected : boolean, found : boolean;
		var u : number, i : number, v : number, j : number;
		var uNode : Nodex, vNode : Nodex;
		var uEdge : Edge, vEdge : Edge;

		undirected = true;

		for(u = 0; u < g.nodes.length; u = u + 1){
			uNode = g.nodes[u];
			for(i = 0; i < uNode.edge.length; i = i + 1){
				uEdge = uNode.edge[i];
				v = uEdge.nodeNr;

				if(u == v){
				}else{
					vNode = g.nodes[v];
					found = false;
					for(j = 0; j < vNode.edge.length && !found; j = j + 1){
						vEdge = vNode.edge[j];

						if(vEdge.nodeNr == u && vEdge.weight == uEdge.weight){
							found = true;
						}
					}

					if(!found){
						undirected = false;
					}
				}
			}
		}

		return undirected;
	}


	export function GetGraphComponents(g : DirectedGraph, componentMembership : NumberArrayReference) : boolean{
		var valid : boolean, found : boolean;
		var i : number, nodeNr : number, componentNr : number, done : number, startNode : number;
		var componentMembershipSet : boolean [];
		var list : NumberArrayReference;

		valid = DirectedGraphIsValid(g) && IsUndirected(g);

		if(valid){
			componentMembership.numberArray = CreateNumberArray(g.nodes.length, 0);
			componentMembershipSet = CreateBooleanArray(g.nodes.length, false);
			list = new NumberArrayReference();
			componentNr = 0;
			done = 0;
			startNode = 0;

			for(; done != g.nodes.length; ){
				/* Find a node not currently in a component.*/
				found = false;
				for(i = 0; i < g.nodes.length && !found; i = i + 1){
					if(!componentMembershipSet[i]){
						startNode = i;
						found = true;
					}
				}

				/* Use DFS to find a component*/
				DepthFirstSearch(g, startNode, list);

				/* Record the component.*/
				for(i = 0; i < list.numberArray.length; i = i + 1){
					nodeNr = list.numberArray[i];
					if(!componentMembershipSet[nodeNr]){
						componentMembership.numberArray[nodeNr] = componentNr;
						componentMembershipSet[nodeNr] = true;
						done = done + 1;
					}
				}

				componentNr = componentNr + 1;
			}

			componentMembershipSet = undefined;
		}

		return valid;
	}


	export function TopologicalSort(g : DirectedGraph, list : NumberArrayReference) : boolean{
		var visited : boolean [];
		var ll : LinkedListNumbers;
		var i : number;
		var valid : boolean;

		valid = !DirectedGraphContainsCycleDFS(g);

		if(valid){

			visited = CreateBooleanArray(g.nodes.length, false);
			ll = CreateLinkedListNumbers();

			for(i = 0; i < g.nodes.length; i = i + 1){
				if(!visited[i]){
					TopologicalSortRecursive(g, g.nodes[i], i, visited, ll);
				}
			}

			list.numberArray = LinkedListNumbersToArray(ll);
			FreeLinkedListNumbers(ll);

			for(i = 0; i < list.numberArray.length/2; i = i + 1){
				SwapElementsOfArray(list.numberArray, i, list.numberArray.length - i - 1);
			}
		}

		return valid;
	}


	export function TopologicalSortRecursive(g : DirectedGraph, node : Nodex, nodeNr : number, visited : boolean [], list : LinkedListNumbers) : void{
		var i : number;
		var e : Edge;

		visited[nodeNr] = true;

		for(i = 0; i < node.edge.length; i = i + 1){
			e = node.edge[i];
			if(!visited[e.nodeNr]){
				TopologicalSortRecursive(g, g.nodes[e.nodeNr], e.nodeNr, visited, list);
			}
		}

		LinkedListAddNumber(list, nodeNr);
	}


	export function ConvertPreviousListToForest(forest : Forest, prev : number []) : void{
		var length : number, i : number, j : number, root : number;
		var found : boolean;

		length = 0;
		for(i = 0; i < prev.length; i = i + 1){
			if(prev[i] == i){
				length = length + 1;
			}
		}

		forest.trees = new Array<Tree>(length);

		j = 0;
		for(i = 0; i < length; i = i + 1){
			/* Find next root.*/
			root = 0;
			found = false;
			for(; j < prev.length && !found; j = j + 1){
				if(prev[j] == j){
					root = j;
					found = true;
				}
			}

			/* Create tree from root.*/
			forest.trees[i] = ConvertPreviousListToTree(root, prev);
		}
	}


	export function ConvertPreviousListToTree(root : number, prev : number []) : Tree{
		var tree : Tree;
		var branches : number, i : number, branch : number;

		tree = new Tree();
		tree.label = root;

		/* Count branches.*/
		branches = 0;
		for(i = 0; i < prev.length; i = i + 1){
			if(prev[i] == root && root != i){
				branches = branches + 1;
			}
		}

		/* Add branches.*/
		tree.branches = new Array<Tree>(branches);
		branch = 0;
		for(i = 0; i < prev.length; i = i + 1){
			if(prev[i] == root && root != i){
				tree.branches[branch] = ConvertPreviousListToTree(i, prev);
				branch = branch + 1;
			}
		}

		return tree;
	}


	export function ConvertLinkedListTreesToForest(forest : Forest, roots : LinkedListNumbers, trees : LinkedListNumbers []) : void{
		var node : LinkedListNodeNumbers;
		var length : number, current : number;

		node = roots.first;
		length = LinkedListNumbersLength(roots);
		forest.trees = new Array<Tree>(length);

		for(current = 0; !node.end; current = current + 1){
			forest.trees[current] = ConvertLinkedListTreeToTree(node.value, trees);
			node = node.next;
		}
	}


	export function ConvertLinkedListTreeToTree(root : number, trees : LinkedListNumbers []) : Tree{
		var tree : Tree;
		var rootList : LinkedListNumbers;
		var node : LinkedListNodeNumbers;
		var current : number, length : number;

		rootList = trees[root];

		tree = new Tree();
		tree.label = root;
		length = LinkedListNumbersLength(rootList);
		tree.branches = new Array<Tree>(length);

		node = rootList.first;

		for(current = 0; !node.end; current = current + 1){
			tree.branches[current] = ConvertLinkedListTreeToTree(node.value, trees);
			node = node.next;
		}

		return tree;
	}


	export function SpanningTreeAlgorithmsTest(failures : NumberReference) : void{
		var g : DirectedGraph;
		var valid : boolean;
		var forest : Forest;

		forest = new Forest();

		g = MakeUndirectedGraphForMST();
		valid = PrimsAlgorithmNoQueue(g, forest);
		CheckUndirectedGraphForMST(failures, valid, forest);
		valid = PrimsAlgorithm(g, forest);
		CheckUndirectedGraphForMST(failures, valid, forest);
		valid = KruskalsAlgorithm(g, forest);
		CheckUndirectedGraphForMSTKruskals(failures, valid, forest);

		g = MakeUndirectedGraph();
		valid = PrimsAlgorithm(g, forest);
		CheckUndirectedGraphMST(failures, valid, forest);
		valid = PrimsAlgorithmNoQueue(g, forest);
		CheckUndirectedGraphMST(failures, valid, forest);
		valid = KruskalsAlgorithm(g, forest);
		CheckUndirectedGraphMSTKruskals(failures, valid, forest);

		g = MakeUndirectedGraphWithThreeComponents();
		valid = PrimsAlgorithm(g, forest);
		CheckUndirectedGraphWithThreeComponentsMST(failures, valid, forest);
		valid = PrimsAlgorithmNoQueue(g, forest);
		CheckUndirectedGraphWithThreeComponentsMST(failures, valid, forest);
		valid = KruskalsAlgorithm(g, forest);
		CheckUndirectedGraphWithThreeComponentsMSTKruskals(failures, valid, forest);

		g = MakeUndirectedGraphForMST2();
		valid = PrimsAlgorithmNoQueue(g, forest);
		/*CheckUndirectedGraphForMST2(failures, valid, forest);*/
		valid = PrimsAlgorithm(g, forest);
		CheckUndirectedGraphForMST2(failures, valid, forest);
		valid = KruskalsAlgorithm(g, forest);
		CheckUndirectedGraphForMST2Kruskals(failures, valid, forest);
	}


	export function CheckUndirectedGraphWithThreeComponentsMST(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 3, failures);
		AssertEquals(forest.trees[0].label, 0, failures);
		AssertEquals(forest.trees[0].branches[0].label, 1, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 2, failures);
		AssertEquals(forest.trees[1].label, 3, failures);
		AssertEquals(forest.trees[2].label, 4, failures);
		AssertEquals(forest.trees[2].branches[0].label, 5, failures);
		AssertEquals(forest.trees[2].branches[0].branches[0].label, 6, failures);
	}


	export function CheckUndirectedGraphWithThreeComponentsMSTKruskals(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 3, failures);
		AssertEquals(forest.trees[0].label, 0, failures);
		AssertEquals(forest.trees[0].branches[0].label, 1, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 2, failures);
		AssertEquals(forest.trees[1].label, 3, failures);
		AssertEquals(forest.trees[2].label, 6, failures);
		AssertEquals(forest.trees[2].branches[0].label, 5, failures);
		AssertEquals(forest.trees[2].branches[0].branches[0].label, 4, failures);
	}


	export function CheckUndirectedGraphMST(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 1, failures);
		AssertEquals(forest.trees[0].label, 0, failures);
		AssertEquals(forest.trees[0].branches[0].label, 1, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 2, failures);
	}


	export function CheckUndirectedGraphMSTKruskals(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 1, failures);
		AssertEquals(forest.trees[0].label, 2, failures);
		AssertEquals(forest.trees[0].branches[0].label, 1, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 0, failures);
	}


	export function CheckUndirectedGraphForMST(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 2, failures);
		AssertEquals(forest.trees[1].label, 4, failures);
		AssertEquals(forest.trees[1].branches.length, 0, failures);
		AssertEquals(forest.trees[0].label, 0, failures);
		AssertEquals(forest.trees[0].branches.length, 2, failures);
		AssertEquals(forest.trees[0].branches[0].label, 3, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 2, failures);
		AssertEquals(forest.trees[0].branches[1].label, 1, failures);
	}


	export function CheckUndirectedGraphForMST2(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 2, failures);
		AssertEquals(forest.trees[1].label, 4, failures);
		AssertEquals(forest.trees[1].branches.length, 0, failures);
		AssertEquals(forest.trees[0].label, 0, failures);
		AssertEquals(forest.trees[0].branches[0].label, 3, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 2, failures);
		AssertEquals(forest.trees[0].branches[1].label, 1, failures);
	}


	export function CheckUndirectedGraphForMSTKruskals(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 2, failures);
		AssertEquals(forest.trees[1].label, 4, failures);
		AssertEquals(forest.trees[1].branches.length, 0, failures);
		AssertEquals(forest.trees[0].label, 0, failures);
		AssertEquals(forest.trees[0].branches.length, 1, failures);
		AssertEquals(forest.trees[0].branches[0].label, 3, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 1, failures);
		AssertEquals(forest.trees[0].branches[0].branches[1].label, 2, failures);
	}


	export function CheckUndirectedGraphForMST2Kruskals(failures : NumberReference, valid : boolean, forest : Forest) : void{
		AssertTrue(valid, failures);
		AssertEquals(forest.trees.length, 2, failures);
		AssertEquals(forest.trees[1].label, 4, failures);
		AssertEquals(forest.trees[1].branches.length, 0, failures);
		AssertEquals(forest.trees[0].label, 2, failures);
		AssertEquals(forest.trees[0].branches[0].label, 3, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].label, 0, failures);
		AssertEquals(forest.trees[0].branches[0].branches[0].branches[0].label, 1, failures);
	}


	export function searchTests(failures : NumberReference) : void{
		DFSTest1(failures);
		DFSTest2(failures);
		DFSTest3(failures);

		BFSTest1(failures);
		BFSTest2(failures);
		BFSTest3(failures);
	}


	export function DFSTest1(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var list : NumberArrayReference;

		directedGraph = MakeGraphWithTwoMixedCycles();

		list = new NumberArrayReference();

		DepthFirstSearch(directedGraph, 0, list);

		AssertEquals(list.numberArray[0], 0, failures);
		AssertEquals(list.numberArray[1], 1, failures);
		AssertEquals(list.numberArray[2], 2, failures);
		AssertEquals(list.numberArray[3], 3, failures);
	}


	export function DFSTest2(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var list : NumberArrayReference;

		directedGraph = MakeGraphForDijkstrasAlgorithm();

		list = new NumberArrayReference();

		DepthFirstSearch(directedGraph, 0, list);

		AssertEquals(list.numberArray[0], 0, failures);
		AssertEquals(list.numberArray[1], 1, failures);
		AssertEquals(list.numberArray[2], 2, failures);
		AssertEquals(list.numberArray[3], 3, failures);
		AssertEquals(list.numberArray[4], 4, failures);
		AssertEquals(list.numberArray[5], 5, failures);
	}


	export function DFSTest3(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var list : NumberArrayReference;

		directedGraph = MakeGraphForDijkstrasAlgorithm2();

		list = CreateNumberArrayReferenceLengthValue(directedGraph.nodes.length, 0);

		DepthFirstSearch(directedGraph, 0, list);

		AssertEquals(list.numberArray[0], 0, failures);
		AssertEquals(list.numberArray[1], 1, failures);
		AssertEquals(list.numberArray[2], 7, failures);
		AssertEquals(list.numberArray[3], 2, failures);
		AssertEquals(list.numberArray[4], 8, failures);
		AssertEquals(list.numberArray[5], 3, failures);
		AssertEquals(list.numberArray[6], 4, failures);
		AssertEquals(list.numberArray[7], 5, failures);
		AssertEquals(list.numberArray[8], 6, failures);
	}


	export function BFSTest1(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var list : NumberArrayReference;

		directedGraph = MakeGraphWithTwoMixedCycles();

		list = new NumberArrayReference();

		BreadthFirstSearch(directedGraph, 0, list);

		AssertEquals(list.numberArray[0], 0, failures);
		AssertEquals(list.numberArray[1], 1, failures);
		AssertEquals(list.numberArray[2], 2, failures);
		AssertEquals(list.numberArray[3], 3, failures);
	}


	export function BFSTest2(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var list : NumberArrayReference;

		directedGraph = MakeGraphForDijkstrasAlgorithm();

		list = new NumberArrayReference();

		BreadthFirstSearch(directedGraph, 0, list);

		AssertEquals(list.numberArray[0], 0, failures);
		AssertEquals(list.numberArray[1], 1, failures);
		AssertEquals(list.numberArray[2], 2, failures);
		AssertEquals(list.numberArray[3], 5, failures);
		AssertEquals(list.numberArray[4], 3, failures);
		AssertEquals(list.numberArray[5], 4, failures);
	}


	export function BFSTest3(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var list : NumberArrayReference;

		directedGraph = MakeGraphForDijkstrasAlgorithm2();

		list = new NumberArrayReference();

		BreadthFirstSearch(directedGraph, 0, list);

		AssertEquals(list.numberArray[0], 0, failures);
		AssertEquals(list.numberArray[1], 1, failures);
		AssertEquals(list.numberArray[2], 7, failures);
		AssertEquals(list.numberArray[3], 6, failures);
		AssertEquals(list.numberArray[4], 2, failures);
		AssertEquals(list.numberArray[5], 8, failures);
		AssertEquals(list.numberArray[6], 5, failures);
		AssertEquals(list.numberArray[7], 3, failures);
		AssertEquals(list.numberArray[8], 4, failures);
	}


	export function ShortestPathsTests(failures : NumberReference) : void{
		testDijkstrasAlgorithm(failures);
		testDijkstrasAlgorithm2(failures);
		testBellmanFordAlgorithm(failures);
		testBellmanFordAlgorithm2(failures);
		testFloydWarshallAlgorithm(failures);
		testFloydWarshallAlgorithm2(failures);
		testShortestPathAlgorithmsDistanceOnly(failures);
	}


	export function testDijkstrasAlgorithm(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var d : NumberArrayReference, p : NumberArrayReference;
		var ds : BooleanArrayReference;

		directedGraph = MakeGraphForDijkstrasAlgorithm();

		d = new NumberArrayReference();
		p = new NumberArrayReference();
		ds = new BooleanArrayReference();
		DijkstrasAlgorithm(directedGraph, 0, d, ds, p);

		CheckShortestPath(failures, d, p);
	}


	export function testShortestPathAlgorithmsDistanceOnly(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var path : NumberArrayReference;
		var distance : NumberReference;
		var success : boolean;

		directedGraph = MakeGraphForDijkstrasAlgorithm();

		path = new NumberArrayReference();
		distance = new NumberReference();
		success = DijkstrasAlgorithmDestinationOnly(directedGraph, 0, 5, path, distance);

		AssertTrue(success, failures);
		AssertEquals(distance.numberValue, 11, failures);
		AssertEquals(path.numberArray.length, 3, failures);
		AssertEquals(path.numberArray[0], 0, failures);
		AssertEquals(path.numberArray[1], 2, failures);
		AssertEquals(path.numberArray[2], 5, failures);

		path = new NumberArrayReference();
		distance = new NumberReference();
		success = BellmanFordAlgorithmDestinationOnly(directedGraph, 0, 5, path, distance);

		AssertTrue(success, failures);
		AssertEquals(distance.numberValue, 11, failures);
		AssertEquals(path.numberArray.length, 3, failures);
		AssertEquals(path.numberArray[0], 0, failures);
		AssertEquals(path.numberArray[1], 2, failures);
		AssertEquals(path.numberArray[2], 5, failures);
	}


	export function testBellmanFordAlgorithm(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var d : NumberArrayReference, p : NumberArrayReference;
		var ds : BooleanArrayReference;
		var success : boolean;

		directedGraph = MakeGraphForDijkstrasAlgorithm();

		d = new NumberArrayReference();
		p = new NumberArrayReference();
		ds = new BooleanArrayReference();
		success = BellmanFordAlgorithm(directedGraph, 0, d, ds, p);

		AssertTrue(success, failures);
		CheckShortestPath(failures, d, p);
	}


	export function CheckShortestPath(failures : NumberReference, d : NumberArrayReference, p : NumberArrayReference) : void{
		AssertEquals(d.numberArray[0], 0, failures);
		AssertEquals(d.numberArray[1], 7, failures);
		AssertEquals(d.numberArray[2], 9, failures);
		AssertEquals(d.numberArray[3], 20, failures);
		AssertEquals(d.numberArray[4], 20, failures);
		AssertEquals(d.numberArray[5], 11, failures);

		AssertEquals(p.numberArray[0], 0, failures);
		AssertEquals(p.numberArray[1], 0, failures);
		AssertEquals(p.numberArray[2], 0, failures);
		AssertEquals(p.numberArray[3], 2, failures);
		AssertEquals(p.numberArray[4], 5, failures);
		AssertEquals(p.numberArray[5], 2, failures);
	}


	export function MakeGraphForDijkstrasAlgorithm() : DirectedGraph{
		var directedGraph : DirectedGraph;

		directedGraph = CreateDirectedGraph(6);

		directedGraph.nodes[0].edge = new Array<Edge>(3);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 7);
		directedGraph.nodes[0].edge[1] = CreateEdge(2, 9);
		directedGraph.nodes[0].edge[2] = CreateEdge(5, 14);

		directedGraph.nodes[1].edge = new Array<Edge>(3);
		directedGraph.nodes[1].edge[0] = CreateEdge(0, 7);
		directedGraph.nodes[1].edge[1] = CreateEdge(2, 10);
		directedGraph.nodes[1].edge[2] = CreateEdge(3, 15);

		directedGraph.nodes[2].edge = new Array<Edge>(4);
		directedGraph.nodes[2].edge[0] = CreateEdge(0, 9);
		directedGraph.nodes[2].edge[1] = CreateEdge(1, 10);
		directedGraph.nodes[2].edge[2] = CreateEdge(3, 11);
		directedGraph.nodes[2].edge[3] = CreateEdge(5, 2);

		directedGraph.nodes[3].edge = new Array<Edge>(3);
		directedGraph.nodes[3].edge[0] = CreateEdge(1, 15);
		directedGraph.nodes[3].edge[1] = CreateEdge(2, 11);
		directedGraph.nodes[3].edge[2] = CreateEdge(4, 6);

		directedGraph.nodes[4].edge = new Array<Edge>(2);
		directedGraph.nodes[4].edge[0] = CreateEdge(3, 6);
		directedGraph.nodes[4].edge[1] = CreateEdge(5, 9);

		directedGraph.nodes[5].edge = new Array<Edge>(3);
		directedGraph.nodes[5].edge[0] = CreateEdge(0, 14);
		directedGraph.nodes[5].edge[1] = CreateEdge(2, 2);
		directedGraph.nodes[5].edge[2] = CreateEdge(4, 9);

		return directedGraph;
	}


	export function testDijkstrasAlgorithm2(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var d : NumberArrayReference, p : NumberArrayReference;
		var ds : BooleanArrayReference;

		directedGraph = MakeGraphForDijkstrasAlgorithm2();
		d = new NumberArrayReference();
		p = new NumberArrayReference();
		ds = new BooleanArrayReference();
		DijkstrasAlgorithm(directedGraph, 0, d, ds, p);

		CheckShortestPath2(failures, d, p);
	}


	export function testBellmanFordAlgorithm2(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var d : NumberArrayReference, p : NumberArrayReference;
		var ds : BooleanArrayReference;

		directedGraph = MakeGraphForDijkstrasAlgorithm2();
		d = new NumberArrayReference();
		p = new NumberArrayReference();
		ds = new BooleanArrayReference();
		BellmanFordAlgorithm(directedGraph, 0, d, ds, p);

		CheckShortestPath2(failures, d, p);
	}


	export function CheckShortestPath2(failures : NumberReference, d : NumberArrayReference, p : NumberArrayReference) : void{
		AssertEquals(d.numberArray[0], 0, failures);
		AssertEquals(d.numberArray[1], 5, failures);
		AssertEquals(d.numberArray[2], 7, failures);
		AssertEquals(d.numberArray[3], 7, failures);
		AssertEquals(d.numberArray[4], 7, failures);
		AssertEquals(d.numberArray[5], 6, failures);
		AssertEquals(d.numberArray[6], 3, failures);
		AssertEquals(d.numberArray[7], 4, failures);
		AssertEquals(d.numberArray[8], 5, failures);

		AssertEquals(p.numberArray[0], 0, failures);
		AssertEquals(p.numberArray[1], 0, failures);
		AssertTrue(p.numberArray[2] == 7 || p.numberArray[2] == 1, failures);
		AssertEquals(p.numberArray[3], 8, failures);
		AssertEquals(p.numberArray[4], 5, failures);
		AssertEquals(p.numberArray[5], 7, failures);
		AssertEquals(p.numberArray[6], 0, failures);
		AssertEquals(p.numberArray[7], 6, failures);
		AssertEquals(p.numberArray[8], 7, failures);
	}


	export function MakeGraphForDijkstrasAlgorithm2() : DirectedGraph{
		var directedGraph : DirectedGraph;

		directedGraph = CreateDirectedGraph(9);

		directedGraph.nodes[0].edge = new Array<Edge>(3);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 5);
		directedGraph.nodes[0].edge[1] = CreateEdge(7, 7);
		directedGraph.nodes[0].edge[2] = CreateEdge(6, 3);

		directedGraph.nodes[1].edge = new Array<Edge>(3);
		directedGraph.nodes[1].edge[0] = CreateEdge(0, 5);
		directedGraph.nodes[1].edge[1] = CreateEdge(7, 3);
		directedGraph.nodes[1].edge[2] = CreateEdge(2, 2);

		directedGraph.nodes[2].edge = new Array<Edge>(4);
		directedGraph.nodes[2].edge[0] = CreateEdge(1, 2);
		directedGraph.nodes[2].edge[1] = CreateEdge(7, 3);
		directedGraph.nodes[2].edge[2] = CreateEdge(8, 3);
		directedGraph.nodes[2].edge[3] = CreateEdge(3, 4);

		directedGraph.nodes[3].edge = new Array<Edge>(3);
		directedGraph.nodes[3].edge[0] = CreateEdge(2, 4);
		directedGraph.nodes[3].edge[1] = CreateEdge(8, 2);
		directedGraph.nodes[3].edge[2] = CreateEdge(4, 5);

		directedGraph.nodes[4].edge = new Array<Edge>(3);
		directedGraph.nodes[4].edge[0] = CreateEdge(3, 5);
		directedGraph.nodes[4].edge[1] = CreateEdge(8, 3);
		directedGraph.nodes[4].edge[2] = CreateEdge(5, 1);

		directedGraph.nodes[5].edge = new Array<Edge>(4);
		directedGraph.nodes[5].edge[0] = CreateEdge(4, 1);
		directedGraph.nodes[5].edge[1] = CreateEdge(8, 2);
		directedGraph.nodes[5].edge[2] = CreateEdge(7, 2);
		directedGraph.nodes[5].edge[3] = CreateEdge(6, 7);

		directedGraph.nodes[6].edge = new Array<Edge>(3);
		directedGraph.nodes[6].edge[0] = CreateEdge(0, 3);
		directedGraph.nodes[6].edge[1] = CreateEdge(7, 1);
		directedGraph.nodes[6].edge[2] = CreateEdge(5, 7);

		directedGraph.nodes[7].edge = new Array<Edge>(6);
		directedGraph.nodes[7].edge[0] = CreateEdge(0, 7);
		directedGraph.nodes[7].edge[1] = CreateEdge(1, 3);
		directedGraph.nodes[7].edge[2] = CreateEdge(2, 3);
		directedGraph.nodes[7].edge[3] = CreateEdge(8, 1);
		directedGraph.nodes[7].edge[4] = CreateEdge(5, 2);
		directedGraph.nodes[7].edge[5] = CreateEdge(6, 1);

		directedGraph.nodes[8].edge = new Array<Edge>(5);
		directedGraph.nodes[8].edge[0] = CreateEdge(3, 2);
		directedGraph.nodes[8].edge[1] = CreateEdge(2, 3);
		directedGraph.nodes[8].edge[2] = CreateEdge(7, 1);
		directedGraph.nodes[8].edge[3] = CreateEdge(5, 2);
		directedGraph.nodes[8].edge[4] = CreateEdge(4, 3);

		return directedGraph;
	}


	export function testFloydWarshallAlgorithm(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var distances : Distances;
		var success : boolean;
		var d : NumberArrayReference, p : NumberArrayReference;
		var i : number;
		var path : number [];

		directedGraph = MakeGraphForDijkstrasAlgorithm();

		distances = CreateDistancesFloydWarshallAlgorithm(directedGraph.nodes.length);
		success = FloydWarshallAlgorithm(directedGraph, distances);

		AssertTrue(success, failures);

		d = CreateNumberArrayReferenceLengthValue(directedGraph.nodes.length, 0);
		p = CreateNumberArrayReferenceLengthValue(directedGraph.nodes.length, 0);

		for(i = 0; i < directedGraph.nodes.length; i = i + 1){
			d.numberArray[i] = distances.fromx[0].to[i].length;
			path = GetPathFromDistances(distances, 0, i);
			if(path.length >= 2){
				p.numberArray[i] = path[path.length - 2];
			}else{
				p.numberArray[i] = i;
			}
		}

		CheckShortestPath(failures, d, p);
	}


	export function testFloydWarshallAlgorithm2(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var distances : Distances;
		var success : boolean;
		var d : NumberArrayReference, p : NumberArrayReference;
		var i : number;
		var path : number [];

		directedGraph = MakeGraphForDijkstrasAlgorithm2();

		distances = CreateDistancesFloydWarshallAlgorithm(directedGraph.nodes.length);
		success = FloydWarshallAlgorithm(directedGraph, distances);

		AssertTrue(success, failures);

		d = CreateNumberArrayReferenceLengthValue(directedGraph.nodes.length, 0);
		p = CreateNumberArrayReferenceLengthValue(directedGraph.nodes.length, 0);

		for(i = 0; i < directedGraph.nodes.length; i = i + 1){
			d.numberArray[i] = distances.fromx[0].to[i].length;
			path = GetPathFromDistances(distances, 0, i);
			if(path.length >= 2){
				p.numberArray[i] = path[path.length - 2];
			}else{
				p.numberArray[i] = i;
			}
		}

		CheckShortestPath2(failures, d, p);
	}


	export function test() : number{
		var failures : NumberReference;

		failures = CreateNumberReference(0);

		testOneCycle(failures);
		testNoCycle(failures);
		testTwoCycles(failures);
		testOneSelfCycle(failures);
		testTwoMixedCycles(failures);
		testMatrixFormConversions(failures);
		ShortestPathsTests(failures);
		searchTests(failures);
		IsUndirectedTests(failures);
		GraphComponentsTest(failures);
		SpanningTreeAlgorithmsTest(failures);
		TestTopologicalSort(failures);

		return failures.numberValue;
	}


	export function testOneSelfCycle(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var valid : boolean, cycle : boolean;
		var cycleCount : number;
		var cycles : Cycle [];

		directedGraph = MakeGraphWithOneSelfCycle();
		valid = DirectedGraphIsValid(directedGraph);
		AssertTrue(valid, failures);
		cycle = DirectedGraphContainsCycleDFS(directedGraph);
		AssertTrue(cycle, failures);
		cycleCount = DirectedGraphCountCyclesDFS(directedGraph);
		AssertEquals(1, cycleCount, failures);
		cycles = DirectedGraphGetCyclesDFS(directedGraph);
		AssertEquals(1, cycles.length, failures);
		AssertEquals(1, cycles[0].nodeNrs.length, failures);
		AssertEquals(2, cycles[0].nodeNrs[0], failures);
	}


	export function testTwoCycles(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var valid : boolean, cycle : boolean;
		var cycleCount : number;
		var cycles : Cycle [];

		directedGraph = MakeGraphWithTwoCycles();
		valid = DirectedGraphIsValid(directedGraph);
		AssertTrue(valid, failures);
		cycle = DirectedGraphContainsCycleDFS(directedGraph);
		AssertTrue(cycle, failures);
		cycleCount = DirectedGraphCountCyclesDFS(directedGraph);
		AssertEquals(2, cycleCount, failures);
		cycles = DirectedGraphGetCyclesDFS(directedGraph);
		AssertEquals(2, cycles.length, failures);
		AssertEquals(3, cycles[0].nodeNrs.length, failures);
		AssertEquals(1, cycles[0].nodeNrs[0], failures);
		AssertEquals(2, cycles[0].nodeNrs[1], failures);
		AssertEquals(3, cycles[0].nodeNrs[2], failures);
		AssertEquals(2, cycles[1].nodeNrs.length, failures);
		AssertEquals(0, cycles[1].nodeNrs[0], failures);
		AssertEquals(4, cycles[1].nodeNrs[1], failures);
	}


	export function testNoCycle(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var valid : boolean, cycle : boolean;
		var cycleCount : number;
		var cycles : Cycle [];

		directedGraph = MakeGraphWithoutCycle();
		valid = DirectedGraphIsValid(directedGraph);
		AssertTrue(valid, failures);
		cycle = DirectedGraphContainsCycleDFS(directedGraph);
		AssertFalse(cycle, failures);
		cycleCount = DirectedGraphCountCyclesDFS(directedGraph);
		AssertEquals(0, cycleCount, failures);
		cycles = DirectedGraphGetCyclesDFS(directedGraph);
		AssertEquals(0, cycles.length, failures);
	}


	export function testOneCycle(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var valid : boolean, cycle : boolean;
		var cycleCount : number;
		var cycles : Cycle [];

		directedGraph = MakeGraphWithOneCycle();
		valid = DirectedGraphIsValid(directedGraph);
		AssertTrue(valid, failures);
		cycle = DirectedGraphContainsCycleDFS(directedGraph);
		AssertTrue(cycle, failures);
		cycleCount = DirectedGraphCountCyclesDFS(directedGraph);
		AssertEquals(1, cycleCount, failures);
		cycles = DirectedGraphGetCyclesDFS(directedGraph);
		AssertEquals(1, cycles.length, failures);
		AssertEquals(3, cycles[0].nodeNrs.length, failures);
		AssertEquals(1, cycles[0].nodeNrs[0], failures);
		AssertEquals(2, cycles[0].nodeNrs[1], failures);
		AssertEquals(3, cycles[0].nodeNrs[2], failures);
	}


	export function testTwoMixedCycles(failures : NumberReference) : void{
		var directedGraph : DirectedGraph;
		var valid : boolean, cycle : boolean;
		var cycleCount : number;
		var cycles : Cycle [];

		directedGraph = MakeGraphWithTwoMixedCycles();
		valid = DirectedGraphIsValid(directedGraph);
		AssertTrue(valid, failures);
		cycle = DirectedGraphContainsCycleDFS(directedGraph);
		AssertTrue(cycle, failures);
		cycleCount = DirectedGraphCountCyclesDFS(directedGraph);
		AssertEquals(2, cycleCount, failures);
		cycles = DirectedGraphGetCyclesDFS(directedGraph);
		AssertEquals(2, cycles.length, failures);
		AssertEquals(3, cycles[0].nodeNrs.length, failures);
		AssertEquals(1, cycles[0].nodeNrs[0], failures);
		AssertEquals(2, cycles[0].nodeNrs[1], failures);
		AssertEquals(3, cycles[0].nodeNrs[2], failures);
		AssertEquals(1, cycles[1].nodeNrs.length, failures);
		AssertEquals(3, cycles[1].nodeNrs[0], failures);
	}


	export function testMatrixFormConversions(failures : NumberReference) : void{
		var a : DirectedGraph, b : DirectedGraph;
		var m1 : DirectedGraphMatrix, m2 : DirectedGraphMatrix;

		a = MakeGraphWithTwoMixedCycles();
		b = CreateDirectedGraphFromMatrixForm(CreateDirectedGraphMatrixFromListForm(MakeGraphWithTwoMixedCycles()));

		AssertTrue(DirectedGraphsEqual(a, b), failures);

		m1 = CreateDirectedGraphMatrixFromListForm(a);
		m2 = CreateDirectedGraphMatrixFromListForm(a);

		AssertTrue(DirectedGraphMatricesEqual(m1, m2), failures);
	}


	export function IsUndirectedTests(failures : NumberReference) : void{
		var g : DirectedGraph;
		var undirected : boolean;

		g = MakeGraphWithTwoMixedCycles();

		undirected = IsUndirected(g);
		AssertFalse(undirected, failures);

		g = MakeGraphWithOneSelfCycle();

		undirected = IsUndirected(g);
		AssertFalse(undirected, failures);

		g = MakeGraphWithOneCycle();

		undirected = IsUndirected(g);
		AssertFalse(undirected, failures);

		g = MakeGraphWithoutCycle();

		undirected = IsUndirected(g);
		AssertFalse(undirected, failures);

		g = MakeGraphWithTwoCycles();

		undirected = IsUndirected(g);
		AssertFalse(undirected, failures);

		g = MakeGraphWithTwoMixedCycles();

		undirected = IsUndirected(g);
		AssertFalse(undirected, failures);

		g = MakeUndirectedGraph();

		undirected = IsUndirected(g);
		AssertTrue(undirected, failures);

		g = MakeUndirectedGraphForMST();

		undirected = IsUndirected(g);
		AssertTrue(undirected, failures);
	}


	export function GraphComponentsTest(failures : NumberReference) : void{
		var g : DirectedGraph;
		var valid : boolean;
		var components : NumberArrayReference;

		g = MakeUndirectedGraph();

		components = new NumberArrayReference();
		valid = GetGraphComponents(g, components);
		AssertTrue(valid, failures);

		AssertEquals(components.numberArray[0], 0, failures);
		AssertEquals(components.numberArray[1], 0, failures);
		AssertEquals(components.numberArray[2], 0, failures);

		g = MakeUndirectedGraphWithThreeComponents();

		components = new NumberArrayReference();
		valid = GetGraphComponents(g, components);
		AssertTrue(valid, failures);

		AssertEquals(components.numberArray[0], 0, failures);
		AssertEquals(components.numberArray[1], 0, failures);
		AssertEquals(components.numberArray[2], 0, failures);
		AssertEquals(components.numberArray[3], 1, failures);
		AssertEquals(components.numberArray[4], 2, failures);
		AssertEquals(components.numberArray[5], 2, failures);
		AssertEquals(components.numberArray[6], 2, failures);

		FreeNumberArrayReference(components);
	}


	export function TestTopologicalSort(failures : NumberReference) : void{
		var g : DirectedGraph;
		var valid : boolean;
		var list : NumberArrayReference;

		g = MakeTopologicalSortGraph();
		list = new NumberArrayReference();

		valid = TopologicalSort(g, list);

		AssertTrue(valid, failures);
		AssertEquals(list.numberArray[0], 5, failures);
		AssertEquals(list.numberArray[1], 4, failures);
		AssertEquals(list.numberArray[2], 2, failures);
		AssertEquals(list.numberArray[3], 3, failures);
		AssertEquals(list.numberArray[4], 1, failures);
		AssertEquals(list.numberArray[5], 0, failures);
	}


	export function MakeGraphWithOneSelfCycle() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(4);

		directedGraph.nodes[0].edge = new Array<Edge>(1);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(2);
		directedGraph.nodes[1].edge[0] = CreateEdge(2, 1);
		directedGraph.nodes[1].edge[1] = CreateEdge(3, 1);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(2, 1);

		directedGraph.nodes[3].edge = new Array<Edge>(0);

		return directedGraph;
	}


	export function MakeGraphWithOneCycle() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(4);

		directedGraph.nodes[0].edge = new Array<Edge>(1);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(1);
		directedGraph.nodes[1].edge[0] = CreateEdge(2, 1);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(3, 1);

		directedGraph.nodes[3].edge = new Array<Edge>(1);
		directedGraph.nodes[3].edge[0] = CreateEdge(1, 1);

		return directedGraph;
	}


	export function MakeGraphWithoutCycle() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(3);

		directedGraph.nodes[0].edge = new Array<Edge>(1);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(1);
		directedGraph.nodes[1].edge[0] = CreateEdge(2, 1);

		directedGraph.nodes[2].edge = new Array<Edge>(0);

		return directedGraph;
	}


	export function MakeGraphWithTwoCycles() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(5);

		directedGraph.nodes[0].edge = new Array<Edge>(2);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 1);
		directedGraph.nodes[0].edge[1] = CreateEdge(4, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(1);
		directedGraph.nodes[1].edge[0] = CreateEdge(2, 1);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(3, 1);

		directedGraph.nodes[3].edge = new Array<Edge>(1);
		directedGraph.nodes[3].edge[0] = CreateEdge(1, 1);

		directedGraph.nodes[4].edge = new Array<Edge>(1);
		directedGraph.nodes[4].edge[0] = CreateEdge(0, 1);

		return directedGraph;
	}


	export function MakeGraphWithTwoMixedCycles() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(4);

		directedGraph.nodes[0].edge = new Array<Edge>(1);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(1);
		directedGraph.nodes[1].edge[0] = CreateEdge(2, 1);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(3, 1);

		directedGraph.nodes[3].edge = new Array<Edge>(2);
		directedGraph.nodes[3].edge[0] = CreateEdge(1, 1);
		directedGraph.nodes[3].edge[1] = CreateEdge(3, 1);

		return directedGraph;
	}


	export function MakeUndirectedGraph() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(3);

		directedGraph.nodes[0].edge = new Array<Edge>(2);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 1);
		directedGraph.nodes[0].edge[1] = CreateEdge(0, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(2);
		directedGraph.nodes[1].edge[0] = CreateEdge(0, 1);
		directedGraph.nodes[1].edge[1] = CreateEdge(2, 1);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(1, 1);

		return directedGraph;
	}


	export function MakeUndirectedGraphWithThreeComponents() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(7);

		directedGraph.nodes[0].edge = new Array<Edge>(2);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 2);
		directedGraph.nodes[0].edge[1] = CreateEdge(0, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(2);
		directedGraph.nodes[1].edge[0] = CreateEdge(0, 2);
		directedGraph.nodes[1].edge[1] = CreateEdge(2, 1);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(1, 1);

		directedGraph.nodes[3].edge = new Array<Edge>(0);

		directedGraph.nodes[4].edge = new Array<Edge>(2);
		directedGraph.nodes[4].edge[0] = CreateEdge(5, 1);
		directedGraph.nodes[4].edge[1] = CreateEdge(4, 1);

		directedGraph.nodes[5].edge = new Array<Edge>(2);
		directedGraph.nodes[5].edge[0] = CreateEdge(4, 1);
		directedGraph.nodes[5].edge[1] = CreateEdge(6, 1);

		directedGraph.nodes[6].edge = new Array<Edge>(1);
		directedGraph.nodes[6].edge[0] = CreateEdge(5, 1);

		return directedGraph;
	}


	export function MakeUndirectedGraphForMST() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(5);

		directedGraph.nodes[0].edge = new Array<Edge>(2);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 2);
		directedGraph.nodes[0].edge[1] = CreateEdge(3, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(2);
		directedGraph.nodes[1].edge[0] = CreateEdge(0, 2);
		directedGraph.nodes[1].edge[1] = CreateEdge(3, 2);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(3, 3);

		directedGraph.nodes[3].edge = new Array<Edge>(3);
		directedGraph.nodes[3].edge[0] = CreateEdge(0, 1);
		directedGraph.nodes[3].edge[1] = CreateEdge(1, 2);
		directedGraph.nodes[3].edge[2] = CreateEdge(2, 3);

		directedGraph.nodes[4].edge = new Array<Edge>(0);

		return directedGraph;
	}


	export function MakeUndirectedGraphForMST2() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(5);

		directedGraph.nodes[0].edge = new Array<Edge>(2);
		directedGraph.nodes[0].edge[0] = CreateEdge(1, 2);
		directedGraph.nodes[0].edge[1] = CreateEdge(3, 1);

		directedGraph.nodes[1].edge = new Array<Edge>(2);
		directedGraph.nodes[1].edge[0] = CreateEdge(0, 2);
		directedGraph.nodes[1].edge[1] = CreateEdge(3, 4);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(3, 3);

		directedGraph.nodes[3].edge = new Array<Edge>(3);
		directedGraph.nodes[3].edge[0] = CreateEdge(0, 1);
		directedGraph.nodes[3].edge[1] = CreateEdge(1, 4);
		directedGraph.nodes[3].edge[2] = CreateEdge(2, 3);

		directedGraph.nodes[4].edge = new Array<Edge>(0);

		return directedGraph;
	}


	export function MakeTopologicalSortGraph() : DirectedGraph{
		var directedGraph : DirectedGraph;
		directedGraph = CreateDirectedGraph(6);

		directedGraph.nodes[0].edge = new Array<Edge>(0);

		directedGraph.nodes[1].edge = new Array<Edge>(0);

		directedGraph.nodes[2].edge = new Array<Edge>(1);
		directedGraph.nodes[2].edge[0] = CreateEdge(3, 1);

		directedGraph.nodes[3].edge = new Array<Edge>(1);
		directedGraph.nodes[3].edge[0] = CreateEdge(1, 1);

		directedGraph.nodes[4].edge = new Array<Edge>(2);
		directedGraph.nodes[4].edge[0] = CreateEdge(0, 1);
		directedGraph.nodes[4].edge[1] = CreateEdge(1, 1);

		directedGraph.nodes[5].edge = new Array<Edge>(2);
		directedGraph.nodes[5].edge[0] = CreateEdge(0, 1);
		directedGraph.nodes[5].edge[1] = CreateEdge(2, 1);

		return directedGraph;
	}


	export function StringToNumberArray(stringx : string []) : number []{
		var i : number;
		var array : number [];

		array = new Array<number>(stringx.length);

		for(i = 0; i < stringx.length; i = i + 1){
			array[i] = stringx[i].charCodeAt(0);
		}
		return array;
	}


	export function NumberArrayToString(array : number []) : string []{
		var i : number;
		var stringx : string [];

		stringx = new Array<string>(array.length);

		for(i = 0; i < array.length; i = i + 1){
			stringx[i] = String.fromCharCode(array[i]);
		}
		return stringx;
	}


	export function NumberArraysEqual(a : number [], b : number []) : boolean{
		var equal : boolean;
		var i : number;

		equal = true;
		if(a.length == b.length){
			for(i = 0; i < a.length && equal; i = i + 1){
				if(a[i] != b[i]){
					equal = false;
				}
			}
		}else{
			equal = false;
		}

		return equal;
	}


	export function BooleanArraysEqual(a : boolean [], b : boolean []) : boolean{
		var equal : boolean;
		var i : number;

		equal = true;
		if(a.length == b.length){
			for(i = 0; i < a.length && equal; i = i + 1){
				if(a[i] != b[i]){
					equal = false;
				}
			}
		}else{
			equal = false;
		}

		return equal;
	}


	export function StringsEqual(a : string [], b : string []) : boolean{
		var equal : boolean;
		var i : number;

		equal = true;
		if(a.length == b.length){
			for(i = 0; i < a.length && equal; i = i + 1){
				if(a[i] != b[i]){
					equal = false;
				}
			}
		}else{
			equal = false;
		}

		return equal;
	}


	export function FillNumberArray(a : number [], value : number) : void{
		var i : number;

		for(i = 0; i < a.length; i = i + 1){
			a[i] = value;
		}
	}


	export function FillString(a : string [], value : string) : void{
		var i : number;

		for(i = 0; i < a.length; i = i + 1){
			a[i] = value;
		}
	}


	export function FillBooleanArray(a : boolean [], value : boolean) : void{
		var i : number;

		for(i = 0; i < a.length; i = i + 1){
			a[i] = value;
		}
	}


	export function FillNumberArrayRange(a : number [], value : number, fromx : number, to : number) : boolean{
		var i : number, length : number;
		var success : boolean;

		if(fromx >= 0 && fromx <= a.length && to >= 0 && to <= a.length && fromx <= to){
			length = to - fromx;
			for(i = 0; i < length; i = i + 1){
				a[fromx + i] = value;
			}

			success = true;
		}else{
			success = false;
		}

		return success;
	}


	export function FillBooleanArrayRange(a : boolean [], value : boolean, fromx : number, to : number) : boolean{
		var i : number, length : number;
		var success : boolean;

		if(fromx >= 0 && fromx <= a.length && to >= 0 && to <= a.length && fromx <= to){
			length = to - fromx;
			for(i = 0; i < length; i = i + 1){
				a[fromx + i] = value;
			}

			success = true;
		}else{
			success = false;
		}

		return success;
	}


	export function FillStringRange(a : string [], value : string, fromx : number, to : number) : boolean{
		var i : number, length : number;
		var success : boolean;

		if(fromx >= 0 && fromx <= a.length && to >= 0 && to <= a.length && fromx <= to){
			length = to - fromx;
			for(i = 0; i < length; i = i + 1){
				a[fromx + i] = value;
			}

			success = true;
		}else{
			success = false;
		}

		return success;
	}


	export function CopyNumberArray(a : number []) : number []{
		var i : number;
		var n : number [];

		n = new Array<number>(a.length);

		for(i = 0; i < a.length; i = i + 1){
			n[i] = a[i];
		}

		return n;
	}


	export function CopyBooleanArray(a : boolean []) : boolean []{
		var i : number;
		var n : boolean [];

		n = new Array<boolean>(a.length);

		for(i = 0; i < a.length; i = i + 1){
			n[i] = a[i];
		}

		return n;
	}


	export function CopyString(a : string []) : string []{
		var i : number;
		var n : string [];

		n = new Array<string>(a.length);

		for(i = 0; i < a.length; i = i + 1){
			n[i] = a[i];
		}

		return n;
	}


	export function CopyNumberArrayRange(a : number [], fromx : number, to : number, copyReference : NumberArrayReference) : boolean{
		var i : number, length : number;
		var n : number [];
		var success : boolean;

		if(fromx >= 0 && fromx <= a.length && to >= 0 && to <= a.length && fromx <= to){
			length = to - fromx;
			n = new Array<number>(length);

			for(i = 0; i < length; i = i + 1){
				n[i] = a[fromx + i];
			}

			copyReference.numberArray = n;
			success = true;
		}else{
			success = false;
		}

		return success;
	}


	export function CopyBooleanArrayRange(a : boolean [], fromx : number, to : number, copyReference : BooleanArrayReference) : boolean{
		var i : number, length : number;
		var n : boolean [];
		var success : boolean;

		if(fromx >= 0 && fromx <= a.length && to >= 0 && to <= a.length && fromx <= to){
			length = to - fromx;
			n = new Array<boolean>(length);

			for(i = 0; i < length; i = i + 1){
				n[i] = a[fromx + i];
			}

			copyReference.booleanArray = n;
			success = true;
		}else{
			success = false;
		}

		return success;
	}


	export function CopyStringRange(a : string [], fromx : number, to : number, copyReference : StringReference) : boolean{
		var i : number, length : number;
		var n : string [];
		var success : boolean;

		if(fromx >= 0 && fromx <= a.length && to >= 0 && to <= a.length && fromx <= to){
			length = to - fromx;
			n = new Array<string>(length);

			for(i = 0; i < length; i = i + 1){
				n[i] = a[fromx + i];
			}

			copyReference.stringx = n;
			success = true;
		}else{
			success = false;
		}

		return success;
	}


	export function IsLastElement(length : number, index : number) : boolean{
		return index + 1 == length;
	}


	export function CreateNumberArray(length : number, value : number) : number []{
		var array : number [];

		array = new Array<number>(length);
		FillNumberArray(array, value);

		return array;
	}


	export function CreateBooleanArray(length : number, value : boolean) : boolean []{
		var array : boolean [];

		array = new Array<boolean>(length);
		FillBooleanArray(array, value);

		return array;
	}


	export function CreateString(length : number, value : string) : string []{
		var array : string [];

		array = new Array<string>(length);
		FillString(array, value);

		return array;
	}


	export function SwapElementsOfArray(A : number [], ai : number, bi : number) : void{
		var tmp : number;

		tmp = A[ai];
		A[ai] = A[bi];
		A[bi] = tmp;
	}


	export function Negate(x : number) : number{
		return -x;
	}


	export function Positive(x : number) : number{
		return +x;
	}


	export function Factorial(x : number) : number{
		var i : number, f : number;
		f = 1;

		for(i = 2; i <= x; i = i + 1){
			f = f*i;
		}

		return f;
	}


	export function Round(x : number) : number{
		return Math.floor(x + 0.5);
	}


	export function BankersRound(x : number) : number{
		var r : number;
		if(Absolute(x - Truncate(x)) == 0.5){
			if(!DivisibleBy(Round(x), 2)){
				r = Round(x) - 1;
			}else{
				r = Round(x);
			}
		}else{
			r = Round(x);
		}

		return r;
	}


	export function Ceil(x : number) : number{
		return Math.ceil(x);
	}


	export function Floor(x : number) : number{
		return Math.floor(x);
	}


	export function Truncate(x : number) : number{
		var t : number;
		if(x >= 0){
			t = Math.floor(x);
		}else{
			t = Math.ceil(x);
		}

		return t;
	}


	export function Absolute(x : number) : number{
		return Math.abs(x);
	}


	export function NaturalLogarithm(x : number) : number{
		return Math.log(x);
	}


	export function Sin(x : number) : number{
		return Math.sin(x);
	}


	export function Cos(x : number) : number{
		return Math.cos(x);
	}


	export function Tan(x : number) : number{
		return Math.tan(x);
	}


	export function Asin(x : number) : number{
		return Math.asin(x);
	}


	export function Acos(x : number) : number{
		return Math.acos(x);
	}


	export function Atan(x : number) : number{
		return Math.atan(x);
	}


	export function Atan2(y : number, x : number) : number{
		var a : number;
		a = 0;

		if(x > 0){
			a = Atan(y/x);
		}else if(x < 0 && y >= 0){
			a = Atan(y/x) + Math.PI;
		}else if(x < 0 && y < 0){
			a = Atan(y/x) - Math.PI;
		}else if(x == 0 && y > 0){
			a = Math.PI/2;
		}else if(x == 0 && y < 0){
			a = -Math.PI/2;
		}

		return a;
	}


	export function Squareroot(x : number) : number{
		return Math.sqrt(x);
	}


	export function Exp(x : number) : number{
		return Math.exp(x);
	}


	export function DivisibleBy(a : number, b : number) : boolean{
		return ((a%b) == 0);
	}


	export function Combinations(n : number, k : number) : number{
		return Factorial(n)/(Factorial(n - k)*Factorial(k));
	}


	export function EpsilonCompare(a : number, b : number, epsilon : number) : boolean{
		return Math.abs(a - b) < epsilon;
	}


	export function GreatestCommonDivisor(a : number, b : number) : number{
		var t : number;
		for(; b != 0; ){
			t = b;
			b = a%b;
			a = t;
		}

		return a;
	}


	export function IsInteger(a : number) : boolean{
		return (a - Math.floor(a)) == 0;
	}


	export function GreatestCommonDivisorWithCheck(a : number, b : number, gcdReference : NumberReference) : boolean{
		var success : boolean;
		var gcd : number;
		if(IsInteger(a) && IsInteger(b)){
			gcd = GreatestCommonDivisor(a, b);
			gcdReference.numberValue = gcd;
			success = true;
		}else{
			success = false;
		}

		return success;
	}


	export function LeastCommonMultiple(a : number, b : number) : number{
		var lcm : number;
		if(a > 0 && b > 0){
			lcm = Math.abs(a*b)/GreatestCommonDivisor(a, b);
		}else{
			lcm = 0;
		}

		return lcm;
	}


	export function Sign(a : number) : number{
		var s : number;
		if(a > 0){
			s = 1;
		}else if(a < 0){
			s = -1;
		}else{
			s = 0;
		}

		return s;
	}


	export function Max(a : number, b : number) : number{
		return Math.max(a, b);
	}


	export function Min(a : number, b : number) : number{
		return Math.min(a, b);
	}


	export function Power(a : number, b : number) : number{
		return a**b;
	}


	export function CreateBooleanReference(value : boolean) : BooleanReference{
		var ref : BooleanReference;
		ref = new BooleanReference();
		ref.booleanValue = value;

		return ref;
	}


	export function CreateBooleanArrayReference(value : boolean []) : BooleanArrayReference{
		var ref : BooleanArrayReference;
		ref = new BooleanArrayReference();
		ref.booleanArray = value;

		return ref;
	}


	export function CreateBooleanArrayReferenceLengthValue(length : number, value : boolean) : BooleanArrayReference{
		var ref : BooleanArrayReference;
		var i : number;
		ref = new BooleanArrayReference();
		ref.booleanArray = new Array<boolean>(length);

		for(i = 0; i < length; i = i + 1){
			ref.booleanArray[i] = value;
		}

		return ref;
	}


	export function FreeBooleanArrayReference(booleanArrayReference : BooleanArrayReference) : void{
		delete booleanArrayReference.booleanArray;
		booleanArrayReference = undefined;
	}


	export function CreateCharacterReference(value : string) : CharacterReference{
		var ref : CharacterReference;
		ref = new CharacterReference();
		ref.characterValue = value;

		return ref;
	}


	export function CreateNumberReference(value : number) : NumberReference{
		var ref : NumberReference;
		ref = new NumberReference();
		ref.numberValue = value;

		return ref;
	}


	export function CreateNumberArrayReference(value : number []) : NumberArrayReference{
		var ref : NumberArrayReference;
		ref = new NumberArrayReference();
		ref.numberArray = value;

		return ref;
	}


	export function CreateNumberArrayReferenceLengthValue(length : number, value : number) : NumberArrayReference{
		var ref : NumberArrayReference;
		var i : number;
		ref = new NumberArrayReference();
		ref.numberArray = new Array<number>(length);

		for(i = 0; i < length; i = i + 1){
			ref.numberArray[i] = value;
		}

		return ref;
	}


	export function FreeNumberArrayReference(numberArrayReference : NumberArrayReference) : void{
		delete numberArrayReference.numberArray;
		numberArrayReference = undefined;
	}


	export function CreateStringReference(value : string []) : StringReference{
		var ref : StringReference;
		ref = new StringReference();
		ref.stringx = value;

		return ref;
	}


	export function CreateStringReferenceLengthValue(length : number, value : string) : StringReference{
		var ref : StringReference;
		var i : number;
		ref = new StringReference();
		ref.stringx = new Array<string>(length);

		for(i = 0; i < length; i = i + 1){
			ref.stringx[i] = value;
		}

		return ref;
	}


	export function FreeStringReference(stringReference : StringReference) : void{
		delete stringReference.stringx;
		stringReference = undefined;
	}


	export function CreateStringArrayReference(strings : StringReference []) : StringArrayReference{
		var ref : StringArrayReference;
		ref = new StringArrayReference();
		ref.stringArray = strings;

		return ref;
	}


	export function CreateStringArrayReferenceLengthValue(length : number, value : string []) : StringArrayReference{
		var ref : StringArrayReference;
		var i : number;
		ref = new StringArrayReference();
		ref.stringArray = new Array<StringReference>(length);

		for(i = 0; i < length; i = i + 1){
			ref.stringArray[i] = CreateStringReference(value);
		}

		return ref;
	}


	export function FreeStringArrayReference(stringArrayReference : StringArrayReference) : void{
		var i : number;
		for(i = 0; i < stringArrayReference.stringArray.length; i = i + 1){
			delete stringArrayReference.stringArray[i];
		}
		delete stringArrayReference.stringArray;
		stringArrayReference = undefined;
	}


	export function CreatePriorityQueueBTNumbers() : PriorityQueueBTNumbers{
		var q : PriorityQueueBTNumbers;

		q = new PriorityQueueBTNumbers();
		q.heap = CreateDynamicArrayNumbers();

		return q;
	}


	export function FreePriorityQueueBTNumbers(q : PriorityQueueBTNumbers) : void{
		FreeDynamicArrayNumbers(q.heap);
		q = undefined;
	}


	export function PeekPriorityQueueBTNumbers(q : PriorityQueueBTNumbers, keyReference : NumberReference) : boolean{
		var found : boolean;

		if(!IsEmptyPriorityQueueBTNumbers(q)){
			keyReference.numberValue = DynamicArrayNumbersIndex(q.heap, 0);
			found = true;
		}else{
			found = false;
		}

		return found;
	}


	export function InsertIntoPriorityQueueBTNumbers(q : PriorityQueueBTNumbers, key : number) : void{
		DynamicArrayAddNumber(q.heap, key);

		if(SizeOfPriorityQueueBTNumbers(q) >= 2){
			SiftUpPriorityQueueBTNumbers(q, q.heap.length - 1);
		}
	}


	export function PopPriorityQueueBTNumbers(q : PriorityQueueBTNumbers, keyReference : NumberReference) : boolean{
		var found : boolean;

		found = PeekPriorityQueueBTNumbers(q, keyReference);

		if(found){
			DeleteTopPriorityQueueBTNumbers(q);
		}

		return found;
	}


	export function DeleteTopPriorityQueueBTNumbers(q : PriorityQueueBTNumbers) : boolean{
		var found : boolean;
		var last : number;

		found = IsEmptyPriorityQueueBTNumbers(q);

		if(!IsEmptyPriorityQueueBTNumbers(q)){
			last = q.heap.length - 1;
			SwapElementsOfArray(q.heap.array, 0, last);

			DynamicArrayRemoveNumber(q.heap, last);

			SiftDownPriorityQueueBTNumbers(q, 0);
		}

		return found;
	}


	export function ArrayToPriorityQueueBTNumbers(keys : number []) : PriorityQueueBTNumbers{
		var q : PriorityQueueBTNumbers;
		var i : number;

		q = CreatePriorityQueueBTNumbers();

		for(i = 0; i < keys.length; i = i + 1){
			InsertIntoPriorityQueueBTNumbers(q, keys[i]);
		}

		return q;
	}


	export function SizeOfPriorityQueueBTNumbers(q : PriorityQueueBTNumbers) : number{
		return q.heap.length;
	}


	export function IsEmptyPriorityQueueBTNumbers(q : PriorityQueueBTNumbers) : boolean{
		return q.heap.length == 0;
	}


	export function SiftUpPriorityQueueBTNumbers(q : PriorityQueueBTNumbers, index : number) : void{
		var parent : number;
		var iKey : number, pKey : number;
		var done : boolean;

		done = false;
		for(; !done && index != 0; ){
			parent = Math.floor((index - 1)/2);

			iKey = DynamicArrayNumbersIndex(q.heap, index);
			pKey = DynamicArrayNumbersIndex(q.heap, parent);

			if(iKey > pKey){
				SwapElementsOfArray(q.heap.array, index, parent);
			}else{
				done = true;
			}

			index = parent;
		}
	}


	export function SiftDownPriorityQueueBTNumbers(q : PriorityQueueBTNumbers, index : number) : void{
		var parent : number, c1 : number, c2 : number;
		var c1Key : number, c2Key : number, pKey : number, size : number;
		var done : boolean;

		size = SizeOfPriorityQueueBTNumbers(q);

		done = false;
		for(; !done; ){
			parent = index;
			c1 = 2*parent + 1;
			c2 = 2*parent + 2;

			pKey = DynamicArrayNumbersIndex(q.heap, parent);
			c1Key = DynamicArrayNumbersIndex(q.heap, c1);
			c2Key = DynamicArrayNumbersIndex(q.heap, c2);

			if(c1Key > pKey && c1 < size || c2Key > pKey && c2 < size){
				if(c1Key >= c2Key && c1 < size){
					SwapElementsOfArray(q.heap.array, c1, parent);
					index = c1;
				}else if(c1Key <= c2Key && c2 < size){
					SwapElementsOfArray(q.heap.array, c2, parent);
					index = c2;
				}else{
					done = true;
				}
			}else{
				done = true;
			}
		}
	}


	export function CreatePriorityQueueBTNumKeyValue() : PriorityQueueBTNumKeyValue{
		var q : PriorityQueueBTNumKeyValue;

		q = new PriorityQueueBTNumKeyValue();
		q.heapKey = CreateDynamicArrayNumbers();
		q.heapValue = CreateDynamicArrayNumbers();

		return q;
	}


	export function FreePriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue) : void{
		FreeDynamicArrayNumbers(q.heapKey);
		FreeDynamicArrayNumbers(q.heapValue);
		q = undefined;
	}


	export function PeekPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue, keyReference : NumberReference, valueReference : NumberReference) : boolean{
		var found : boolean;

		if(!IsEmptyPriorityQueueBTNumKeyValue(q)){
			keyReference.numberValue = DynamicArrayNumbersIndex(q.heapKey, 0);
			valueReference.numberValue = DynamicArrayNumbersIndex(q.heapValue, 0);
			found = true;
		}else{
			found = false;
		}

		return found;
	}


	export function InsertIntoPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue, key : number, value : number) : void{
		DynamicArrayAddNumber(q.heapKey, key);
		DynamicArrayAddNumber(q.heapValue, value);

		if(SizeOfPriorityQueueBTNumKeyValue(q) >= 2){
			SiftUpPriorityQueueBTNumKeyValue(q, q.heapKey.length - 1);
		}
	}


	export function PopPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue, keyReference : NumberReference, valueReference : NumberReference) : boolean{
		var found : boolean;

		found = PeekPriorityQueueBTNumKeyValue(q, keyReference, valueReference);

		if(found){
			DeleteTopPriorityQueueBTNumKeyValue(q);
		}

		return found;
	}


	export function DeleteTopPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue) : boolean{
		var found : boolean;
		var last : number;

		found = IsEmptyPriorityQueueBTNumKeyValue(q);

		if(!IsEmptyPriorityQueueBTNumKeyValue(q)){
			last = q.heapKey.length - 1;
			SwapElementsOfArray(q.heapKey.array, 0, last);
			SwapElementsOfArray(q.heapValue.array, 0, last);

			DynamicArrayRemoveNumber(q.heapKey, last);
			DynamicArrayRemoveNumber(q.heapValue, last);

			SiftDownPriorityQueueBTNumKeyValue(q, 0);
		}

		return found;
	}


	export function ArrayToPriorityQueueBTNumKeyValue(keys : number [], values : number []) : PriorityQueueBTNumKeyValue{
		var q : PriorityQueueBTNumKeyValue;
		var i : number;

		q = CreatePriorityQueueBTNumKeyValue();

		for(i = 0; i < keys.length; i = i + 1){
			InsertIntoPriorityQueueBTNumKeyValue(q, keys[i], values[i]);
		}

		return q;
	}


	export function SizeOfPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue) : number{
		return q.heapKey.length;
	}


	export function IsEmptyPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue) : boolean{
		return q.heapKey.length == 0;
	}


	export function SiftUpPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue, index : number) : void{
		var parent : number;
		var iKey : number, pKey : number;
		var done : boolean;

		done = false;
		for(; !done && index != 0; ){
			parent = Math.floor((index - 1)/2);

			iKey = DynamicArrayNumbersIndex(q.heapKey, index);
			pKey = DynamicArrayNumbersIndex(q.heapKey, parent);

			if(iKey > pKey){
				SwapElementsOfArray(q.heapKey.array, index, parent);
				SwapElementsOfArray(q.heapValue.array, index, parent);
			}else{
				done = true;
			}

			index = parent;
		}
	}


	export function SiftDownPriorityQueueBTNumKeyValue(q : PriorityQueueBTNumKeyValue, index : number) : void{
		var parent : number, c1 : number, c2 : number;
		var c1Key : number, c2Key : number, pKey : number, size : number;
		var done : boolean;

		size = SizeOfPriorityQueueBTNumKeyValue(q);

		c1Key = 0;
		c2Key = 0;
		done = false;
		for(; !done; ){
			parent = index;
			c1 = 2*parent + 1;
			c2 = 2*parent + 2;

			pKey = DynamicArrayNumbersIndex(q.heapKey, parent);
			if(c1 < size){
				c1Key = DynamicArrayNumbersIndex(q.heapKey, c1);
			}
			if(c2 < size){
				c2Key = DynamicArrayNumbersIndex(q.heapKey, c2);
			}

			if(c1Key > pKey && c1 < size || c2Key > pKey && c2 < size){
				if(c1Key >= c2Key && c1 < size){
					SwapElementsOfArray(q.heapKey.array, c1, parent);
					SwapElementsOfArray(q.heapValue.array, c1, parent);
					index = c1;
				}else if(c1Key <= c2Key && c2 < size){
					SwapElementsOfArray(q.heapKey.array, c2, parent);
					SwapElementsOfArray(q.heapValue.array, c2, parent);
					index = c2;
				}else{
					done = true;
				}
			}else{
				done = true;
			}
		}
	}


	export function CreateLinkedListNumbers() : LinkedListNumbers{
		var ll : LinkedListNumbers;

		ll = new LinkedListNumbers();
		ll.first = new LinkedListNodeNumbers();
		ll.last = ll.first;
		ll.last.end = true;

		return ll;
	}


	export function CreateLinkedListNumbersArray(length : number) : LinkedListNumbers []{
		var lls : LinkedListNumbers [];
		var i : number;

		lls = new Array<LinkedListNumbers>(length);
		for(i = 0; i < lls.length; i = i + 1){
			lls[i] = CreateLinkedListNumbers();
		}

		return lls;
	}


	export function LinkedListAddNumber(ll : LinkedListNumbers, value : number) : void{
		ll.last.end = false;
		ll.last.value = value;
		ll.last.next = new LinkedListNodeNumbers();
		ll.last.next.end = true;
		ll.last = ll.last.next;
	}


	export function LinkedListNumbersLength(ll : LinkedListNumbers) : number{
		var l : number;
		var node : LinkedListNodeNumbers;

		l = 0;
		node = ll.first;
		for(; !node.end; ){
			node = node.next;
			l = l + 1;
		}

		return l;
	}


	export function LinkedListNumbersIndex(ll : LinkedListNumbers, index : number) : number{
		var i : number;
		var node : LinkedListNodeNumbers;

		node = ll.first;
		for(i = 0; i < index; i = i + 1){
			node = node.next;
		}

		return node.value;
	}


	export function LinkedListInsertNumber(ll : LinkedListNumbers, index : number, value : number) : void{
		var i : number;
		var node : LinkedListNodeNumbers, tmp : LinkedListNodeNumbers;

		if(index == 0){
			tmp = ll.first;
			ll.first = new LinkedListNodeNumbers();
			ll.first.next = tmp;
			ll.first.value = value;
			ll.first.end = false;
		}else{
			node = ll.first;
			for(i = 0; i < index - 1; i = i + 1){
				node = node.next;
			}

			tmp = node.next;
			node.next = new LinkedListNodeNumbers();
			node.next.next = tmp;
			node.next.value = value;
			node.next.end = false;
		}
	}


	export function LinkedListSet(ll : LinkedListNumbers, index : number, value : number) : void{
		var i : number;
		var node : LinkedListNodeNumbers;

		node = ll.first;
		for(i = 0; i < index; i = i + 1){
			node = node.next;
		}

		node.next.value = value;
	}


	export function LinkedListRemoveNumber(ll : LinkedListNumbers, index : number) : void{
		var i : number;
		var node : LinkedListNodeNumbers, prev : LinkedListNodeNumbers;

		node = ll.first;
		prev = ll.first;

		for(i = 0; i < index; i = i + 1){
			prev = node;
			node = node.next;
		}

		if(index == 0){
			ll.first = prev.next;
		}

		prev.next = prev.next.next;
	}


	export function FreeLinkedListNumbers(ll : LinkedListNumbers) : void{
		var node : LinkedListNodeNumbers, prev : LinkedListNodeNumbers;

		node = ll.first;

		for(; !node.end; ){
			prev = node;
			node = node.next;
			prev = undefined;
		}

		node = undefined;
	}


	export function FreeLinkedListNumbersArray(lls : LinkedListNumbers []) : void{
		var i : number;

		for(i = 0; i < lls.length; i = i + 1){
			FreeLinkedListNumbers(lls[i]);
		}
		lls = undefined;
	}


	export function LinkedListNumbersToArray(ll : LinkedListNumbers) : number []{
		var array : number [];
		var length : number, i : number;
		var node : LinkedListNodeNumbers;

		node = ll.first;

		length = LinkedListNumbersLength(ll);

		array = new Array<number>(length);

		for(i = 0; i < length; i = i + 1){
			array[i] = node.value;
			node = node.next;
		}

		return array;
	}


	export function ArrayToLinkedListNumbers(array : number []) : LinkedListNumbers{
		var ll : LinkedListNumbers;
		var i : number;

		ll = CreateLinkedListNumbers();

		for(i = 0; i < array.length; i = i + 1){
			LinkedListAddNumber(ll, array[i]);
		}

		return ll;
	}


	export function LinkedListNumbersEqual(a : LinkedListNumbers, b : LinkedListNumbers) : boolean{
		var equal : boolean, done : boolean;
		var an : LinkedListNodeNumbers, bn : LinkedListNodeNumbers;

		an = a.first;
		bn = b.first;

		equal = true;
		done = false;
		for(; equal && !done; ){
			if(an.end == bn.end){
				if(an.end){
					done = true;
				}else if(an.value == bn.value){
					an = an.next;
					bn = bn.next;
				}else{
					equal = false;
				}
			}else{
				equal = false;
			}
		}

		return equal;
	}


	export function CreateDynamicArrayNumbers() : DynamicArrayNumbers{
		var da : DynamicArrayNumbers;

		da = new DynamicArrayNumbers();
		da.array = new Array<number>(10);
		da.length = 0;

		return da;
	}


	export function CreateDynamicArrayNumbersWithInitialCapacity(capacity : number) : DynamicArrayNumbers{
		var da : DynamicArrayNumbers;

		da = new DynamicArrayNumbers();
		da.array = new Array<number>(capacity);
		da.length = 0;

		return da;
	}


	export function DynamicArrayAddNumber(da : DynamicArrayNumbers, value : number) : void{
		if(da.length == da.array.length){
			DynamicArrayNumbersIncreaseSize(da);
		}

		da.array[da.length] = value;
		da.length = da.length + 1;
	}


	export function DynamicArrayNumbersIncreaseSize(da : DynamicArrayNumbers) : void{
		var newLength : number, i : number;
		var newArray : number [];

		newLength = Math.round(da.array.length*3/2);
		newArray = new Array<number>(newLength);

		for(i = 0; i < da.array.length; i = i + 1){
			newArray[i] = da.array[i];
		}

		delete da.array;

		da.array = newArray;
	}


	export function DynamicArrayNumbersDecreaseSizeNecessary(da : DynamicArrayNumbers) : boolean{
		var needsDecrease : boolean;

		needsDecrease = false;

		if(da.length > 10){
			needsDecrease = da.length <= Math.round(da.array.length*2/3);
		}

		return needsDecrease;
	}


	export function DynamicArrayNumbersDecreaseSize(da : DynamicArrayNumbers) : void{
		var newLength : number, i : number;
		var newArray : number [];

		newLength = Math.round(da.array.length*2/3);
		newArray = new Array<number>(newLength);

		for(i = 0; i < da.array.length; i = i + 1){
			newArray[i] = da.array[i];
		}

		delete da.array;

		da.array = newArray;
	}


	export function DynamicArrayNumbersIndex(da : DynamicArrayNumbers, index : number) : number{
		return da.array[index];
	}


	export function DynamicArrayNumbersLength(da : DynamicArrayNumbers) : number{
		return da.length;
	}


	export function DynamicArrayInsertNumber(da : DynamicArrayNumbers, index : number, value : number) : void{
		var i : number;

		if(da.length == da.array.length){
			DynamicArrayNumbersIncreaseSize(da);
		}

		for(i = da.length; i >= index; i = i - 1){
			da.array[i + 1] = da.array[i];
		}

		da.array[index] = value;

		da.length = da.length + 1;
	}


	export function DynamicArraySet(da : DynamicArrayNumbers, index : number, value : number) : void{
		da.array[index] = value;
	}


	export function DynamicArrayRemoveNumber(da : DynamicArrayNumbers, index : number) : void{
		var i : number;

		for(i = index; i < da.length - 1; i = i + 1){
			da.array[i] = da.array[i + 1];
		}

		da.length = da.length - 1;

		if(DynamicArrayNumbersDecreaseSizeNecessary(da)){
			DynamicArrayNumbersDecreaseSize(da);
		}
	}


	export function FreeDynamicArrayNumbers(da : DynamicArrayNumbers) : void{
		delete da.array;
		da = undefined;
	}


	export function DynamicArrayNumbersToArray(da : DynamicArrayNumbers) : number []{
		var array : number [];
		var i : number;

		array = new Array<number>(da.length);

		for(i = 0; i < da.length; i = i + 1){
			array[i] = da.array[i];
		}

		return array;
	}


	export function ArrayToDynamicArrayNumbersWithOptimalSize(array : number []) : DynamicArrayNumbers{
		var da : DynamicArrayNumbers;
		var i : number;
		var c : number, n : number, newCapacity : number;

		/*
         c = 10*(3/2)^n
         log(c) = log(10*(3/2)^n)
         log(c) = log(10) + log((3/2)^n)
         log(c) = 1 + log((3/2)^n)
         log(c) - 1 = log((3/2)^n)
         log(c) - 1 = n*log(3/2)
         n = (log(c) - 1)/log(3/2)
        */
		c = array.length;
		n = (Math.log(c) - 1)/Math.log(3/2);
		newCapacity = Math.floor(n) + 1;

		da = CreateDynamicArrayNumbersWithInitialCapacity(newCapacity);

		for(i = 0; i < array.length; i = i + 1){
			da.array[i] = array[i];
		}

		return da;
	}


	export function ArrayToDynamicArrayNumbers(array : number []) : DynamicArrayNumbers{
		var da : DynamicArrayNumbers;

		da = new DynamicArrayNumbers();
		da.array = CopyNumberArray(array);
		da.length = array.length;

		return da;
	}


	export function DynamicArrayNumbersEqual(a : DynamicArrayNumbers, b : DynamicArrayNumbers) : boolean{
		var equal : boolean;
		var i : number;

		equal = true;
		if(a.length == b.length){
			for(i = 0; i < a.length && equal; i = i + 1){
				if(a.array[i] != b.array[i]){
					equal = false;
				}
			}
		}else{
			equal = false;
		}

		return equal;
	}


	export function DynamicArrayNumbersToLinkedList(da : DynamicArrayNumbers) : LinkedListNumbers{
		var ll : LinkedListNumbers;
		var i : number;

		ll = CreateLinkedListNumbers();

		for(i = 0; i < da.length; i = i + 1){
			LinkedListAddNumber(ll, da.array[i]);
		}

		return ll;
	}


	export function LinkedListToDynamicArrayNumbers(ll : LinkedListNumbers) : DynamicArrayNumbers{
		var da : DynamicArrayNumbers;
		var i : number;
		var node : LinkedListNodeNumbers;

		node = ll.first;

		da = new DynamicArrayNumbers();
		da.length = LinkedListNumbersLength(ll);

		da.array = new Array<number>(da.length);

		for(i = 0; i < da.length; i = i + 1){
			da.array[i] = node.value;
			node = node.next;
		}

		return da;
	}


	export function AddNumber(list : number [], a : number) : number []{
		var newlist : number [];
		var i : number;

		newlist = new Array<number>(list.length + 1);
		for(i = 0; i < list.length; i = i + 1){
			newlist[i] = list[i];
		}
		newlist[list.length] = a;

		list = undefined;

		return newlist;
	}


	export function AddNumberRef(list : NumberArrayReference, i : number) : void{
		list.numberArray = AddNumber(list.numberArray, i);
	}


	export function RemoveNumber(list : number [], n : number) : number []{
		var newlist : number [];
		var i : number;

		newlist = new Array<number>(list.length - 1);

		if(n >= 0 && n < list.length){
			for(i = 0; i < list.length; i = i + 1){
				if(i < n){
					newlist[i] = list[i];
				}
				if(i > n){
					newlist[i - 1] = list[i];
				}
			}

			list = undefined;
		}else{
			newlist = undefined;
		}

		return newlist;
	}


	export function GetNumberRef(list : NumberArrayReference, i : number) : number{
		return list.numberArray[i];
	}


	export function RemoveNumberRef(list : NumberArrayReference, i : number) : void{
		list.numberArray = RemoveNumber(list.numberArray, i);
	}


	export function AddString(list : StringReference [], a : StringReference) : StringReference []{
		var newlist : StringReference [];
		var i : number;

		newlist = new Array<StringReference>(list.length + 1);

		for(i = 0; i < list.length; i = i + 1){
			newlist[i] = list[i];
		}
		newlist[list.length] = a;

		list = undefined;

		return newlist;
	}


	export function AddStringRef(list : StringArrayReference, i : StringReference) : void{
		list.stringArray = AddString(list.stringArray, i);
	}


	export function RemoveString(list : StringReference [], n : number) : StringReference []{
		var newlist : StringReference [];
		var i : number;

		newlist = new Array<StringReference>(list.length - 1);

		if(n >= 0 && n < list.length){
			for(i = 0; i < list.length; i = i + 1){
				if(i < n){
					newlist[i] = list[i];
				}
				if(i > n){
					newlist[i - 1] = list[i];
				}
			}

			list = undefined;
		}else{
			newlist = undefined;
		}

		return newlist;
	}


	export function GetStringRef(list : StringArrayReference, i : number) : StringReference{
		return list.stringArray[i];
	}


	export function RemoveStringRef(list : StringArrayReference, i : number) : void{
		list.stringArray = RemoveString(list.stringArray, i);
	}


	export function AddBoolean(list : boolean [], a : boolean) : boolean []{
		var newlist : boolean [];
		var i : number;

		newlist = new Array<boolean>(list.length + 1);
		for(i = 0; i < list.length; i = i + 1){
			newlist[i] = list[i];
		}
		newlist[list.length] = a;

		list = undefined;

		return newlist;
	}


	export function AddBooleanRef(list : BooleanArrayReference, i : boolean) : void{
		list.booleanArray = AddBoolean(list.booleanArray, i);
	}


	export function RemoveBoolean(list : boolean [], n : number) : boolean []{
		var newlist : boolean [];
		var i : number;

		newlist = new Array<boolean>(list.length - 1);

		if(n >= 0 && n < list.length){
			for(i = 0; i < list.length; i = i + 1){
				if(i < n){
					newlist[i] = list[i];
				}
				if(i > n){
					newlist[i - 1] = list[i];
				}
			}

			list = undefined;
		}else{
			newlist = undefined;
		}

		return newlist;
	}


	export function GetBooleanRef(list : BooleanArrayReference, i : number) : boolean{
		return list.booleanArray[i];
	}


	export function RemoveDecimalRef(list : BooleanArrayReference, i : number) : void{
		list.booleanArray = RemoveBoolean(list.booleanArray, i);
	}


	export function AddCharacter(list : string [], a : string) : string []{
		var newlist : string [];
		var i : number;

		newlist = new Array<string>(list.length + 1);
		for(i = 0; i < list.length; i = i + 1){
			newlist[i] = list[i];
		}
		newlist[list.length] = a;

		list = undefined;

		return newlist;
	}


	export function AddCharacterRef(list : StringReference, i : string) : void{
		list.stringx = AddCharacter(list.stringx, i);
	}


	export function RemoveCharacter(list : string [], n : number) : string []{
		var newlist : string [];
		var i : number;

		newlist = new Array<string>(list.length - 1);

		if(n >= 0 && n < list.length){
			for(i = 0; i < list.length; i = i + 1){
				if(i < n){
					newlist[i] = list[i];
				}
				if(i > n){
					newlist[i - 1] = list[i];
				}
			}

			list = undefined;
		}else{
			newlist = undefined;
		}

		return newlist;
	}


	export function GetCharacterRef(list : StringReference, i : number) : string{
		return list.stringx[i];
	}


	export function RemoveCharacterRef(list : StringReference, i : number) : void{
		list.stringx = RemoveCharacter(list.stringx, i);
	}


	export function TreeHeight(tree : Tree) : number{
		var height : number, i : number, branchHeight : number;
		var heightSet : boolean;

		heightSet = false;
		height = 0;

		for(i = 0; i < tree.branches.length; i = i + 1){
			branchHeight = TreeHeight(tree.branches[i]);
			if(!heightSet){
				height = branchHeight;
				heightSet = true;
			}else if(branchHeight > height){
				height = branchHeight;
			}
		}

		if(tree.branches.length == 0){
			height = 0;
		}else{
			height = height + 1;
		}

		return height;
	}


	export function TreeNumberOfNodes(tree : Tree) : number{
		var nodes : number, i : number;

		nodes = 0;

		for(i = 0; i < tree.branches.length; i = i + 1){
			nodes = nodes + TreeNumberOfNodes(tree.branches[i]);
		}

		return nodes + 1;
	}


	export function AssertFalse(b : boolean, failures : NumberReference) : void{
		if(b){
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertTrue(b : boolean, failures : NumberReference) : void{
		if(!b){
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertEquals(a : number, b : number, failures : NumberReference) : void{
		if(a != b){
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertBooleansEqual(a : boolean, b : boolean, failures : NumberReference) : void{
		if(a != b){
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertCharactersEqual(a : string, b : string, failures : NumberReference) : void{
		if(a != b){
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertStringEquals(a : string [], b : string [], failures : NumberReference) : void{
		if(!StringsEqual(a, b)){
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertNumberArraysEqual(a : number [], b : number [], failures : NumberReference) : void{
		var i : number;

		if(a.length == b.length){
			for(i = 0; i < a.length; i = i + 1){
				AssertEquals(a[i], b[i], failures);
			}
		}else{
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertBooleanArraysEqual(a : boolean [], b : boolean [], failures : NumberReference) : void{
		var i : number;

		if(a.length == b.length){
			for(i = 0; i < a.length; i = i + 1){
				AssertBooleansEqual(a[i], b[i], failures);
			}
		}else{
			failures.numberValue = failures.numberValue + 1;
		}
	}


	export function AssertStringArraysEqual(a : StringReference [], b : StringReference [], failures : NumberReference) : void{
		var i : number;

		if(a.length == b.length){
			for(i = 0; i < a.length; i = i + 1){
				AssertStringEquals(a[i].stringx, b[i].stringx, failures);
			}
		}else{
			failures.numberValue = failures.numberValue + 1;
		}
	}


