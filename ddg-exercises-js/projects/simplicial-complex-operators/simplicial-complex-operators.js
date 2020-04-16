"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                mesh.vertices.forEach( (vert, idx) => {vert.index = idx;})
                mesh.edges.forEach( (edge, idx) => {edge.index = idx;})
                mesh.faces.forEach( (face, idx) => {face.index = idx;})  
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                let A = new Triplet(mesh.edges.length, mesh.vertices.length);

                mesh.edges.forEach( ( edge )=> {
                        const {vertex: v1, twin: {vertex: v2}} = edge.halfedge;
                        A.addEntry(1, edge.index, v1.index);
                        A.addEntry(1, edge.index, v2.index);
                })

                return SparseMatrix.fromTriplet(A);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let A = new Triplet(mesh.faces.length, mesh.edges.length);
                mesh.faces.forEach( ( face )=> {
                         const {edge: e1, next: {edge: e2, next: {edge: e3}}} = face.halfedge;
                         A.addEntry(1, face.index, e1.index);
                         A.addEntry(1, face.index, e2.index);
                         A.addEntry(1, face.index, e3.index);
                })
                 return SparseMatrix.fromTriplet(A);
              
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let res = DenseMatrix.zeros(this.mesh.vertices.length, 1);
                
                subset.vertices.forEach((idx) => { 
                        res.set(1,idx,0);
                });

                return res;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let res = DenseMatrix.zeros(this.mesh.edges.length, 1);
                subset.edges.forEach((idx) => { 
                        res.set(1,idx,0);
                });

                return res;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let res = DenseMatrix.zeros(this.mesh.faces.length, 1);
                subset.faces.forEach((idx) => { 
                        res.set(1,idx,0);
                });

                return res;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                let res = MeshSubset.deepCopy(subset),
                    subsetVertices = this.buildVertexVector(subset);
                    
        
                
                for(var i = 0; i < subsetVertices.nRows(); i++){
                        if (subsetVertices.get(i,0) == 1){
                            var edges = new VertexEdgeIterator(this.mesh.vertices[i].halfedge);
                            for(let edge of edges){
                                res.addEdge(edge.index);
                            }
                           
                        }
                };
                var faceVect = this.A1.timesDense(this.buildEdgeVector(res));
                for(var r=0;r<faceVect.nRows();r++)
                        if(faceVect.get(r,0)>0) res.addFace(r);

                return res;
               
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                let res = MeshSubset.deepCopy(subset),
                subsetEdges = this.buildEdgeVector(subset),
                subsetFaces = this.buildFaceVector(subset);
                
                for(var idx = 0; idx < subsetEdges.nRows(); idx++){
                        if (subsetEdges.get(idx,0) == 1){
                         const { vertex: {index: v0},  next: {vertex: {index: v1}} } = this.mesh.edges[idx].halfedge;
                         res.addVertex(v0);
                         res.addVertex(v1);
                        }
                };
        
        
                for(var idx = 0; idx < subsetFaces.nRows(); idx++){
                        if (subsetFaces.get(idx,0) == 1){
                         let vertices = this.mesh.faces[idx].adjacentVertices(),
                             edges = this.mesh.faces[idx].adjacentEdges();

                                for (let v of vertices) {
                                        res.addVertex(v.index);
                        
                                }
                                for (let e of edges) {
                                        res.addEdge(e.index);
                        
                                }       
                        }
                }

                return res; 
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                let comp = this.closure(this.star(subset)),
                    invComp = this.star(this.closure(subset));
                
                comp.deleteVertices(invComp.vertices);
                comp.deleteEdges(invComp.edges);
                comp.deleteFaces(invComp.faces);

                return comp; 
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                return (this.closure(subset).equals(subset));
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                var degree = 0;

                let facesSize = subset.faces.size,
                    subsetEdges = this.buildEdgeVector(subset),
                    faceEdgeVector = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                    
                if (facesSize > 0 ){
                        for(var i=0;i<subsetEdges.nRows();i++)
				{  
                                        var se = subsetEdges.get(i,0);
					var fe = faceEdgeVector.get(i,0);
					if((fe != se) && (se || fe)){
						return -1;
                                        }
                        }
                        degree = 2;
                }

                let     edgesSize = subset.edges.size,
                        subsetVertices = this.buildVertexVector(subset),
                        edgeVertVector = this.A0.transpose().timesDense(this.buildEdgeVector(subset));

                if (edgesSize > 0){
                        for(var i=0;i<subsetVertices.nRows();i++)
				{  
                                        var sv = subsetVertices.get(i,0);
                                        var ev = edgeVertVector.get(i,0);
                                        //expetec to be equivalent (ev != sv) && (sv|| ev)
					if(((sv!=0)!=(ev!=0))){
						return -1;
                                        }
                        }
                       if (degree!= 2) degree = 1;
                }

                return degree;


        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */

        boundary(subset) {
                var res = new MeshSubset()

                var faceEdgeVector = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                for(var i=0;i<faceEdgeVector.nRows();i++)
                {
                        if(faceEdgeVector.get(i,0)==1) res.addEdge(i);
                }
                
                if(res.edges.size==0) 
			{
                        var edgeVertVector = this.A0.transpose().timesDense(this.buildEdgeVector(subset));
                        for(var i=0;i<edgeVertVector.nRows();i++)
                        {
                                if(edgeVertVector.get(i,0)==1) res.addVertex(i);
                        }
                        
		}else{
                        for (let eidx of res.edges) {
                                const { vertex: {index: v0},  next: {vertex: {index: v1}} } = this.mesh.edges[eidx].halfedge;
                                 res.addVertex(v0);
                                 res.addVertex(v1);
                        }
        

                }

                return res;
        }
}
