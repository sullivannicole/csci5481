from skbio.tree import TreeNode as tn

####
# Example 1: build example tree pictured in HW3 Newick explanation, root down
####

print("\nBuilding example tree from the root down...")

# start with root, for example
tr = tn('F') # root node
tr.children += [tn('A', .1, parent=tr)] # add child to initially empty list
tr.children += [tn('B', .2, parent=tr)] # remember to set parent 
tr.children += [tn('E', .5, parent=tr)]

tr.children[2].children += [tn('C',.3)] # we know "E" is the 3rd child of "F"
tr.children[2].children += [tn('D',.4)]

# We could have used "find" but it would be slower
# tr.find('E').children += [tn('C',.3)]
# tr.find('E').children += [tn('D',.4)]

print(tr)
print(tr.ascii_art())
print("\nPre-order traversal:", end=" ")
for n in tr.preorder():
	print(n.name, end=" ")
print()
print("Post-order traversal:", end=" ")
for n in tr.postorder():
	print(n.name, end=" ")
print()


####
# Example 2: build example tree from HW, tips up
####

print("\nBuilding example tree from the tips up...")

# create all tips
tip1 = tn('A') # tip node
tip2 = tn('B') # tip node
tip3 = tn('C') # tip node
tip4 = tn('D') # tip node

# connect tips "C" and "D" to first internal node "E"
inode1 = tn('E',children=[tip3, tip4]) # internal node E
tip3.parent = inode1 # don't forget to set parent nodes
tip4.parent = inode1
tip3.length = .3
tip4.length = .4

# connect tips "A" and "B" and internal node "E" to internal node "F"
inode2 = tn('F',children=[tip1, tip2, inode1]) # internal node "F"
tip1.parent = inode2 # don't forget to set parent nodes
tip2.parent = inode2
inode1.parent = inode2
tip1.length = .1
tip2.length = .2
inode1.length = .5

tr = inode2 # use the root node as the tree
print(tr)


# write Newick format file
print("\nWriting Newick format output file tree.tre...")
tr.write('tree.tre')
# write example tip labels file if you want to use the hw3 R script to plot
open('tip-labels.txt','w').write('A\t1\tblue\nB\t2\tgreen\nC\t3\tred\nD\t4\tred\n')
