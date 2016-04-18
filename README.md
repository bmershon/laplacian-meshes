# laplacian-meshes

Please view this README rendered by GitHub at https://github.com/bmershon/shape-google

*All images, words, and code contained in this repository may be reproduced so long as the original author is given credit (Chris Tralie and Brooks Mershon).*

This assignment was completed as part of a course in 3D Digital Geometry (Math 290) taken at Duke University during Spring 2016. The course was taught by [Chris Tralie](http://www.ctralie.com/).

## Laplacian Mesh Editing

### Laplacian Matrix

The Laplacian operator is encoded as a sparse matrix **L**, with anchor rows appended to encode the weights of the anchor vertices (which may be manually moved, hence the name Laplacian *editing*).

### Cotangent Weights

Rather than using equal weights for each neighboring vertex in the Laplacian operator, we can attempt to correct for irregular mesh resolutions by using [Cotangent Weights](http://www.ctralie.com/Teaching/COMPSCI290/Assignments/Group3_LaplacianMesh/).

*Homer's arms are raised by placing two anchors at his finger tips and moving them vertically. A small handful of anchors placed symmetrically about his body help to restrict edits to his arms. Cotangent weighting is used here*

<img src="img/homer.png" width="404">
<img src="img/homer-arms-raised.png" width="404">


