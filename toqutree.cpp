
/**
 *
 * toqutree (pa3)
 * significant modification of a quadtree .
 * toqutree.cpp
 * This file will be used for grading.
 *
 */

#include "toqutree.h"
#include "cmath"

toqutree::Node::Node(pair<int,int> ctr, int dim, HSLAPixel a)
	:center(ctr),dimension(dim),avg(a),NW(NULL),NE(NULL),SE(NULL),SW(NULL)
	{}

toqutree::~toqutree(){
	clear(root);
}

//copy toqutree
toqutree::toqutree(const toqutree & other) {
	root = copy(other.root);
}


toqutree & toqutree::operator=(const toqutree & rhs){
	if (this != &rhs) {
		clear(root);
		root = copy(rhs.root);
	}
	return *this;
}

/* called by destructor and assignment operator */
void toqutree::clear(Node * & curr){
	if (curr == NULL){
		return;
	}
	clear(curr->NW);
	clear(curr->NE);
	clear(curr->SE);
	clear(curr->SW);
	delete(curr);
	curr = NULL;
}

/* called by assignment operator and copy constructor */
toqutree::Node * toqutree::copy(const Node * other) {
	if (other == NULL){
		return NULL;
	}
	Node* newRoot = new Node(other->center, other->dimension, other->avg);
	newRoot->NW = copy(other->NW);
	newRoot->NE = copy(other->NE);
	newRoot->SE = copy(other->SE);
	newRoot->SW = copy(other->SW);
	return newRoot;
}

int toqutree::size() {
	return size(root);
}

int toqutree::size(const Node* root){
	if (root == NULL){
		return 0;
	} else {
		return 1 + size(root->NW) + size(root->NE) + size(root->SE) + size(root->SW);
	}
}

/* This constructor grabs the 2^k x 2^k sub-image centered */
/* in imIn and uses it to build a quadtree. It may assume  */
/* that imIn is large enough to contain an image of that size. */
toqutree::toqutree(PNG & imIn, int k){
	int newDims = pow(2,k);
	int centerX = (imIn.width()/2)-1;
	int centerY = (imIn.height()/2)-1;
	int startX = (imIn.width() - newDims)/2;
	int startY = (imIn.height() - newDims)/2;
	PNG* croppedImg = new PNG(newDims,newDims);
	HSLAPixel* oldPixel, *newPixel;
	if (k == 0){
		HSLAPixel* oldPixel = imIn.getPixel(centerX,centerY);
		root = new Node(make_pair(centerX,centerY), 0, *oldPixel);
	} else{
	for (int i = 0; i < newDims; i++){
		for (int j = 0; j < newDims; j++){
		oldPixel = imIn.getPixel(j+startY, i+startX);
		newPixel = croppedImg->getPixel(j, i);
		*newPixel = *oldPixel;
		}
	}
	root = buildTree(croppedImg, k);
	delete(croppedImg);
	}
}

// Note that you will want to practice careful memory use
// In this function. We pass the dynamically allocated image
// via pointer so that it may be released after it is used .
// similarly, at each level of the tree you will want to
// declare a dynamically allocated stats object, and free it
// once you've used it to choose a split point, and calculate
// an average.
toqutree::Node * toqutree::buildTree(PNG * im, int k) {
	stats* rootStat = new stats(*im);
	HSLAPixel avg = rootStat->getAvg(make_pair(0,0),make_pair(im->width()-1,im->width()-1));

	if (k == 1){
		Node* output = new Node(make_pair(1,1),1,avg);
		output->NW = new Node(make_pair(0,0),1,*im->getPixel(0,0));
		output->NE = new Node(make_pair(0,0),1,*im->getPixel(1,0));
		output->SW = new Node(make_pair(0,0),1,*im->getPixel(0,1));
		output->SE = new Node(make_pair(0,0),1,*im->getPixel(1,1));
		delete(rootStat);
		return output;
	} else {
		int imDim = im->width(); 									//2^k by 2^k, same dimensions
		int smallDim = pow(2,k-1);								//small image dimensions
		int startX = (imDim - smallDim)/2;
		int startY = (imDim - smallDim)/2;
		PNG* NWChild, *NEChild, *SEChild, *SWChild;  //Declare the children
		PNG* quadImage = new PNG(imDim*2,imDim*2);
		makeFourQuad(quadImage, *im, imDim);
		stats quadStats = stats(*quadImage);
		double minEntropy = 9007199254740992;
		pair<int,int> tempCtr(0,0);
		for (int i = startX; i < smallDim+startX; i++){
			for (int j = startY; j < smallDim+startY; j++){
				pair<int,int> SE(j,i);
				pair<int,int> SElr(j+smallDim-1,i+smallDim-1);

				pair<int,int> SW(j+smallDim,i);
				pair<int,int> SWlr(j+smallDim+smallDim-1,i+smallDim-1);

				pair<int,int> NE(j,i+smallDim);
				pair<int,int> NElr(j+smallDim-1,i+smallDim+smallDim-1);

				pair<int,int> NW(j+smallDim,i+smallDim);
				pair<int,int> NWlr(j+smallDim+smallDim-1,i+smallDim+smallDim-1);

				double entropySE = quadStats.entropy(SE,SElr);
				double entropySW = quadStats.entropy(SW,SWlr);
				double entropyNE = quadStats.entropy(NE,NElr);
				double entropyNW = quadStats.entropy(NW,NWlr);
				double avgEntropy = (entropySE+entropySW+entropyNE+entropyNW)/4;
				if (avgEntropy < minEntropy){
					minEntropy = avgEntropy;
					tempCtr.first = SE.first;
					tempCtr.second = SE.second;
				}
			}
		}
		NWChild = new PNG(smallDim,smallDim);
		NEChild = new PNG(smallDim,smallDim);
		SEChild = new PNG(smallDim,smallDim);
		SWChild = new PNG(smallDim,smallDim);
		int ctrX = tempCtr.first;
		int ctrY = tempCtr.second;
		for (int i = 0; i < smallDim; i++){
			for (int j = 0; j < smallDim; j++){
				HSLAPixel* newPix = SEChild->getPixel(j,i);
				HSLAPixel* oldPix = quadImage->getPixel(j+ctrX,i+ctrY);
				*newPix = *oldPix;
				newPix = SWChild->getPixel(j,i);
				oldPix = quadImage->getPixel(j+smallDim+ctrX,i+ctrY);
				*newPix = *oldPix;
				newPix = NEChild->getPixel(j,i);
				oldPix = quadImage->getPixel(j+ctrX,i+smallDim+ctrY);
				*newPix = *oldPix;
				newPix = NWChild->getPixel(j,i);
				oldPix = quadImage->getPixel(j+smallDim+ctrX,i+smallDim+ctrY);
				*newPix = *oldPix;
			}
		}
		Node* output = new Node(tempCtr,k,avg);
		output->NW = buildTree(NWChild,k-1);
		output->NE = buildTree(NEChild,k-1);
		output->SE = buildTree(SEChild,k-1);
		output->SW = buildTree(SWChild,k-1);
		delete(quadImage);
		delete(NWChild);
		delete(NEChild);
		delete(SEChild);
		delete(SWChild);
		delete(rootStat);
		return output;
	}
}

// My algorithm for this problem included a helper function
// that was analogous to Find in a BST, but it navigated the
// quadtree, instead.
PNG toqutree::render(){
	int dim = pow(2, root->dimension);
	PNG output = PNG(dim,dim);
	render(&output, root, root->dimension);
	return output;
	// PNG* temp = new PNG(dim,dim);
	// render(temp, root, root->dimension);
	// PNG output = *temp;
	// delete(temp);
	// return output;
}

/* Render helper to deal with inputs */
void toqutree::render(PNG* img, Node* subRoot, int k){
	int dim = pow(2,k);
  int smallDim = pow(2,k-1);
  if (hasNoChild(subRoot)){
    for (int i = 0; i < dim; i++){
      for (int j = 0; j < dim; j++){
        *img->getPixel(j,i) = subRoot->avg;
      }
    }
  }else {
		int ctrX = subRoot->center.first;
		int ctrY = subRoot->center.second;
    PNG* NWChild = new PNG(smallDim,smallDim);
    PNG* NEChild = new PNG(smallDim,smallDim);
    PNG* SEChild = new PNG(smallDim,smallDim);
    PNG* SWChild = new PNG(smallDim,smallDim);
    render(NWChild, subRoot->NW, k-1);
    render(NEChild, subRoot->NE, k-1);
    render(SEChild, subRoot->SE, k-1);
    render(SWChild, subRoot->SW, k-1);
    for (int i = 0; i < smallDim; i++){
      for (int j = 0; j < smallDim; j++){
        *img->getPixel(modInt(j+ctrX, dim),modInt(i+ctrY, dim))
				= *SEChild->getPixel(j,i);
      }
    }
    for (int i = 0; i < smallDim; i++){
      for (int j = 0; j < smallDim; j++){
        *img->getPixel(modInt(j+ctrX-smallDim, dim),modInt(i+ctrY, dim))
				= *SWChild->getPixel(j,i);
      }
    }
    for (int i = 0; i < smallDim; i++){
      for (int j = 0; j < smallDim; j++){
        *img->getPixel(modInt(j+ctrX, dim),modInt(i+ctrY-smallDim, dim))
				= *NEChild->getPixel(j,i);
      }
    }
    for (int i = 0; i < smallDim; i++){
      for (int j = 0; j < smallDim; j++){
        *img->getPixel(modInt(j+ctrX-smallDim, dim),modInt(i+ctrY-smallDim, dim))
				= *NWChild->getPixel(j,i);
      }
    }
		delete(NWChild);
		delete(NEChild);
		delete(SEChild);
		delete(SWChild);
	}
}

		/* Alternate implementation using a quadImage */

		// PNG* quadImage = new PNG(2*dim,2*dim);
		// makeFourQuad(quadImage, *img, dim);
		// int ctrX = subRoot->center.first;
		// int ctrY = subRoot->center.second;
		// PNG* SEChild = new PNG(smallDim,smallDim);
		// render(SEChild, subRoot->SE, k-1);
		// for (int i = 0; i < smallDim; i++){
		// 	for (int j = 0; j < smallDim; j++){
		// 		*quadImage->getPixel(j+ctrX,i+ctrY) = *SEChild->getPixel(j,i);
		// 	}
		// }
		// delete(SEChild);
		// PNG* SWChild = new PNG(smallDim, smallDim);
		// render(SWChild, subRoot->SW, k-1);
		// for (int i = 0; i < smallDim; i++){
		// 	for (int j = 0; j < smallDim; j++){
		// 		*quadImage->getPixel(j+smallDim+ctrX,i+ctrY) = *SWChild->getPixel(j,i);
		// 	}
		// }
		// delete(SWChild);
		// PNG* NEChild = new PNG(smallDim, smallDim);
		// render(NEChild, subRoot->NE, k-1);
		// for (int i = 0; i < smallDim; i++){
		// 	for (int j = 0; j < smallDim; j++){
		// 		*quadImage->getPixel(j+ctrX,i+smallDim+ctrY) = *NEChild->getPixel(j,i);
		// 	}
		// }
		// delete(NEChild);
		// PNG* NWChild = new PNG(smallDim, smallDim);
		// render(NWChild, subRoot->NW, k-1);
		// for (int i = 0; i < smallDim; i++){
		// 	for (int j = 0; j < smallDim; j++){
		// 		*quadImage->getPixel(j+smallDim+ctrX,i+smallDim+ctrY) = *NWChild->getPixel(j,i);
		// 	}
		// }
		// delete(NWChild);
		// for (int i = 0; i < dim; i++){
		// 	for (int j = 0; j < dim; j++){
		// 		*img->getPixel((j+ctrX) % dim, (i+ctrY) % dim) = *quadImage->getPixel(j+ctrX,i+ctrY);
		// 	}
		// }
		// delete(quadImage);


/* oops, i left the implementation of this one in the file! */
void toqutree::prune(double tol){
	if (root == NULL){
		return;
	}
	prune(root,tol);
}

// memory test passes even if you just set children to NULL
void toqutree::prune(Node* subRoot, double tol){
	if (hasNoChild(subRoot)){
		return;
	}
	if (shouldPrune(subRoot, tol, subRoot->avg)){
		clear(subRoot->NW);
		subRoot->NW = NULL;
		clear(subRoot->NE);
		subRoot->NE = NULL;
		clear(subRoot->SW);
		subRoot->SW = NULL;
		clear(subRoot->SE);
		subRoot->SE = NULL;
	}else {
		prune(subRoot->NW, tol);
		prune(subRoot->NE, tol);
		prune(subRoot->SW, tol);
		prune(subRoot->SE, tol);
	}
}

bool toqutree::shouldPrune(Node* subRoot, double tol, HSLAPixel rootAvg){
	if (hasNoChild(subRoot)){
		return (rootAvg.dist(subRoot->avg) <= tol);
	}
	return (shouldPrune(subRoot->NW,tol,rootAvg) &&
	shouldPrune(subRoot->NE,tol,rootAvg) &&
	shouldPrune(subRoot->SW,tol,rootAvg) &&
	shouldPrune(subRoot->SE,tol,rootAvg));
}

/* checks if this Node has children */
bool toqutree::hasNoChild(Node* curr){
	if ((curr->NW == NULL) && (curr->NE == NULL) && (curr->SE == NULL) && (curr->SW == NULL)){
		return true;
	} else {
		return false;
	}
}

/* Private helper that makes 3 additional copies of image im
 * to place in 4 quadrants, simplifies split calcuation */
void toqutree::makeFourQuad(PNG* output, PNG im, int imSize){
	for (int i = 0; i < imSize*2; i++){
		for (int j = 0; j < imSize*2; j++){
			HSLAPixel* smallPix = im.getPixel(j % imSize, i % imSize);
			HSLAPixel* bigPix = output->getPixel(j,i);
			*bigPix = *smallPix;
		}
	}
}

int toqutree::modInt(int a, int b){
  if (a >= 0){
    return a % b;
  } else{
    return b + a;
  }
}
