package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.ArrayList;

public class SearchAlgorithms {
	// Attribute, if the object is found, success is true
	boolean success;
	// Attribute records all the nodes and paths that go through during the search
	ArrayList<Node> searchResultSequence;
	// Attribute records the whole cost of a certain type of search
	int totalCost = 0;
	// Attribute, the name of of a certain type of search
	String sequence = "";
	// Attribute records the whole cost of the iterative deepening search
	int totalCostForIterate = 0;

	/**
	 * The Breadth First Search
	 */
	public ArrayList<Node> BreadthFirstSearch(DigitalBinaryTree binaryTree,
			int objective) {
		success = false;
		totalCost = 0;
		sequence = "";
		int i = 0;
		searchResultSequence = new ArrayList<Node>();
		// Add the root of the tree
		searchResultSequence.add(binaryTree.getBinaryTree().get(0));
		while (i < searchResultSequence.size()) {
			if (searchResultSequence.get(i).getValue() == objective) {
				success = true;
				break;
			} else {
				if (null != searchResultSequence.get(i).getLeftChild()) {
					searchResultSequence.add(searchResultSequence.get(i).getLeftChild());
					totalCost += searchResultSequence.get(i).getToLeftChildCost();
					if (searchResultSequence.get(i).getLeftChild().getValue() == objective) {
						success = true;
						break;
					}
				}
				if (null != searchResultSequence.get(i).getRightChild()) {
					searchResultSequence.add(searchResultSequence.get(i).getRightChild());
					totalCost += searchResultSequence.get(i).getToRightChildCost();
					if (searchResultSequence.get(i).getRightChild().getValue() == objective) {
						success = true;
						break;
					}
				}
			}
			i++;
		}
		this.searchResultDisplay(this.success, "Breadth_First_Search");
		return searchResultSequence;
	}

	/**
	 * The Depth First Search
	 */
	public ArrayList<Node> DepthFirstSearch(DigitalBinaryTree binaryTree,
			int objective, int methodType) {
		success = false;
		totalCost = 0;
		sequence = "";
		// Data structure acted as a stack that saves the potential nodes to be checked
		ArrayList<Node> searchBufferSequence = new ArrayList<Node>();
		int i = 0;
		int currentNodeIndex = 0;
		searchResultSequence = new ArrayList<Node>();
		// Add the root of the tree
		searchResultSequence.add(binaryTree.getBinaryTree().get(0));
		do {
			if (searchResultSequence.get(i).getValue() == objective) {
				success = true;
				break;
			} else {
				currentNodeIndex = searchResultSequence.get(i).getValue() - 1;
				// If a parent node has both child nodes, then the left child node should be examine next time, yet the right child node
				// should be stored in the potential stack
				if (binaryTree.getBinaryTree().get(currentNodeIndex).getLeftChild() != (null)) {
					searchResultSequence.add(binaryTree.getBinaryTree().get(currentNodeIndex).getLeftChild());
					totalCost += binaryTree.getBinaryTree().get(currentNodeIndex).getToLeftChildCost();
					if (binaryTree.getBinaryTree().get(currentNodeIndex).getLeftChild().getValue() == objective) {
						success = true;
						break;
					}
					if (binaryTree.getBinaryTree().get(currentNodeIndex).getRightChild() != (null)) {
						searchBufferSequence.add(binaryTree.getBinaryTree().get(currentNodeIndex).getRightChild());
					}
				}
				// If a parent node has no child node, then the node in potential stack should be used
				else if (searchBufferSequence.size() > 0) {
					Node readyToVerify = searchBufferSequence.get(searchBufferSequence.size() - 1);
					searchResultSequence.add(readyToVerify);
					totalCost += binaryTree.getBinaryTree().get((int) readyToVerify.getValue() / 2 - 1).getToRightChildCost();
					searchBufferSequence.remove(searchBufferSequence.size() - 1);
					if (searchBufferSequence.size() == 0 && null != readyToVerify.getLeftChild()) {
						searchBufferSequence.add(readyToVerify.getLeftChild());
					}
				}
			}
			i++;
		} while (searchBufferSequence.size() > 0);
		// As to the Iterative Depth Limited Search, its total cost should contain the cost generated in each iteration
		if (methodType == 2) {
			this.totalCostForIterate += totalCost;
			System.out.println("totalCost of this loop -----: " + totalCost);
		}
		if (methodType == 0)
			this.searchResultDisplay(this.success, "Depth_First_Search");
		else if (methodType == 1)
			this.searchResultDisplay(this.success, "Depth_Limited_Search");
		else if (methodType == 2)
			this.searchResultDisplay(this.success, "Iterative_Depth_Limited_Search");

		return searchResultSequence;
	}

	/**
	 * The Depth Limited First Search
	 */
	public ArrayList<Node> DepthLimitedSearch(DigitalBinaryTree binaryTree,
			int objective, int depthLimit) {
		// if the provided depth is larger than the depth of the original tree, then the search process is just the same as depth first search
		if (binaryTree.treeDepth() <= depthLimit) {
			searchResultSequence = this.DepthFirstSearch(binaryTree, objective,1);
		} else {
			// Trim the tree as the given depth, then execute the depth first search on the new tree
			DigitalBinaryTree vitualBinaryTree = binaryTree.treeTrim(binaryTree, depthLimit);
			searchResultSequence = this.DepthFirstSearch(vitualBinaryTree, objective, 1);
		}
		return searchResultSequence;
	}

	/**
	 * The Iterative Depth Limited Search
	 */
	public ArrayList<Node> IterativeDepthLimitedSearch(DigitalBinaryTree binaryTree,
			int objective) {
		// Attribute iterative time, increasing from 1 to the depth of the binaryTree
		int Iterative_times = 1;
		System.out.println("");
		System.out.println("Iterative Depth Limited Search Begins: ");
		while (Iterative_times <= binaryTree.treeDepth()) {
			System.out.println("This is the Iteration " + Iterative_times + ":");
            // Trim the tree as the iterative time, then execute the depth first search on the new tree until the objective is found or all nodes
			// have been checked
			DigitalBinaryTree vitualBinaryTree = binaryTree.treeTrim(binaryTree, Iterative_times);
			searchResultSequence = this.DepthFirstSearch(vitualBinaryTree, objective, 2);
			if (this.success == true) {
				break;
			}
			Iterative_times++;
		}
		return searchResultSequence;
	}

	
	/**
	 * The Uniform Cost Search
	 */
	public ArrayList<Node> UniformCostSearch(DigitalBinaryTree binaryTree,
			int objective) {
		success = false;
		totalCost = 0;
		sequence = "";
		Node loweatCostNode;
		// Data structure acted as a heap that saves the potential nodes to be checked,  each time choose the one has the lowest cost
		ArrayList<Node> searchBufferSequence = new ArrayList<Node>();
		int i = 0;
		searchResultSequence = new ArrayList<Node>();
		// Add the root of the tree
		searchResultSequence.add(binaryTree.getBinaryTree().get(0));
		do {
			if (searchResultSequence.get(i).getValue() == objective) {
				success = true;
				break;
			} else {
					// If a parent node has child(or children), put it(or them) into heap. Then the one which has lowest cost should be examine next time
					if (null != searchResultSequence.get(i).getLeftChild()) {
						searchBufferSequence.add(searchResultSequence.get(i).getLeftChild());	
						if (null != searchResultSequence.get(i).getRightChild()) {				
							searchBufferSequence.add(searchResultSequence.get(i).getRightChild());
						}
					}
					loweatCostNode=new Node();
					loweatCostNode.setTotalReachCost(Integer.MAX_VALUE);
					for(int m=0;m<searchBufferSequence.size();m++){	
						if(searchBufferSequence.get(m).getTotalReachCost()<loweatCostNode.getTotalReachCost()){
							loweatCostNode=searchBufferSequence.get(m);
						}
					}
					searchResultSequence.add(loweatCostNode);
					// The single path cost = TotalReachCost of a node - TotalReachCost of its parent node
					totalCost += (loweatCostNode.getTotalReachCost()-binaryTree.getBinaryTree().get((int)(loweatCostNode.getValue()/2 - 1))
							.getTotalReachCost());
					if (loweatCostNode.getValue() == objective) {
						success = true;
						break;
					}
					searchBufferSequence.remove(loweatCostNode);
			}
			i++;
		} while (searchBufferSequence.size() > 0);
		this.searchResultDisplay(this.success, "Uniform_Cost_Search");

		return searchResultSequence;
	}
	
	/**
	 * A Variant of Uniform Cost Search
	 */
	public ArrayList<Node> VariantUniformCostSearch(DigitalBinaryTree binaryTree,
			int objective) {
		success = false;
		totalCost = 0;
		sequence = "";
		// Data structure acted as a stack that saves the potential nodes to be checked
		ArrayList<Node> searchBufferSequence = new ArrayList<Node>();
		int i = 0;
		searchResultSequence = new ArrayList<Node>();
		// Add the root of the tree
		searchResultSequence.add(binaryTree.getBinaryTree().get(0));
		do {
			if (searchResultSequence.get(i).getValue() == objective) {
				success = true;
				break;
			} else {
				// If a parent node has both child nodes, then the one which have a lower cost should be examine next time, yet the other 
				//should be stored in the potential stack
				if (null != searchResultSequence.get(i).getLeftChild() && null != searchResultSequence.get(i).getRightChild()) {
					if (searchResultSequence.get(i).getToLeftChildCost() <= searchResultSequence.get(i).getToRightChildCost()) {
						searchResultSequence.add(searchResultSequence.get(i).getLeftChild());
						totalCost += searchResultSequence.get(i).getToLeftChildCost();
						if (searchResultSequence.get(i).getLeftChild().getValue() == objective) {
							success = true;
							break;
						}
						searchBufferSequence.add(searchResultSequence.get(i).getRightChild());
					} else {
						searchResultSequence.add(searchResultSequence.get(i).getRightChild());
						totalCost += searchResultSequence.get(i).getToRightChildCost();
						if (searchResultSequence.get(i).getRightChild().getValue() == objective) {
							success = true;
							break;
						}
						searchBufferSequence.add(searchResultSequence.get(i).getLeftChild());
					}
				} else if (null != searchResultSequence.get(i).getLeftChild()) {
					searchResultSequence.add(searchResultSequence.get(i).getLeftChild());
					totalCost += searchResultSequence.get(i).getToLeftChildCost();
					if (searchResultSequence.get(i).getLeftChild().getValue() == objective) {
						success = true;
						break;
					}
				}
				// If a parent node has no child node, then the node in potential stack should be used
				else {
					Node readyToVerify = searchBufferSequence.get(searchBufferSequence.size() - 1);
					searchResultSequence.add(readyToVerify);
					// determine the current node is whether the left child or right child of its parent
					if (readyToVerify.getValue() % 2 == 0) {
						totalCost += binaryTree.getBinaryTree().get((int) readyToVerify.getValue() / 2 - 1).getToLeftChildCost();
					} else {
						totalCost += binaryTree.getBinaryTree().get((int) readyToVerify.getValue() / 2 - 1).getToRightChildCost();
					}
					searchBufferSequence.remove(searchBufferSequence.size() - 1);
					// determine the next node should be put in the potential stack is whether the left child or right child of the current node
					if (searchBufferSequence.size() == 0) {
						if (null != readyToVerify.getLeftChild() && null != readyToVerify.getLeftChild()) {
							if (readyToVerify.getToLeftChildCost() <= readyToVerify.getToRightChildCost())
								searchBufferSequence.add(readyToVerify.getLeftChild());
							else
								searchBufferSequence.add(readyToVerify.getRightChild());
						} else if (null != readyToVerify.getLeftChild())
							searchBufferSequence.add(readyToVerify.getLeftChild());
					}
				}
			}
			i++;
		} while (searchBufferSequence.size() > 0);
		this.searchResultDisplay(this.success, "Variant_Uniform_Cost_Search");

		return searchResultSequence;
	}

	/**
	 * Output all the information of a certain Search
	 */
	public void searchResultDisplay(boolean success, String searchName) {
		if (success) {
			System.out.println();
			if (searchName == "Iterative_Depth_Limited_Search") {
				System.out.println("OK, objective found. The total cost of the "+ searchName + " is ["+ this.totalCostForIterate
						+ "]. Here is the search sequence:");
			} else {
				System.out.println("OK, objective found. The total cost of the "+ searchName + " is [" + totalCost+ "]. Here is the search sequence:");
			}
			for (int m = 0; m < searchResultSequence.size() - 1; m++) {
				sequence += "Node(value="+ searchResultSequence.get(m).getValue() + ") "+ "---> ";
			}
			sequence += "Node(value="+ searchResultSequence.get(searchResultSequence.size() - 1).getValue() + ") ";
			System.out.println(sequence);
		} else {
			System.out.println();
			System.out.println("Sorry, such objective does not exist in the current tree.");
		}
	}
	
	
	public static void main(String args[]) {
		// establish the binary tree, total 31 nodes
		DigitalBinaryTree binaryTree = new DigitalBinaryTree(31);
		binaryTree.display();

		// Invoke different tree search algorithms, with object 30
		SearchAlgorithms searchAlgo = new SearchAlgorithms();
		searchAlgo.BreadthFirstSearch(binaryTree, 30);
		searchAlgo.DepthFirstSearch(binaryTree, 30, 0);
		searchAlgo.DepthLimitedSearch(binaryTree, 30, 5);
		searchAlgo.IterativeDepthLimitedSearch(binaryTree, 30);
		searchAlgo.UniformCostSearch(binaryTree, 30);
		// searchAlgo.VariantUniformCostSearch(binaryTree, 30);
	}
}
