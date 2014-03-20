using System;
using System.Collections.Generic;

namespace Delaunay_Voronoi_Library
{
    class Edge
    {
       
        #region _fields

        private List<Vertex> vertices;

        #endregion _fields
 
        #region _constructors
        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Edge with the specified vertices.
        /// </summary>
        /// <param name="vertex1">The vertex value.</param>
        /// <param name="vertex2">The vertex value.</param>
        public Edge(Vertex vertex1, Vertex vertex2)
        {
            vertices = new List<Vertex>();
            vertices.Add(vertex1);
            vertices.Add(vertex2);
        }
        
        #endregion _constructors

        #region _properties
        /// <summary>
        /// Gets the segment vertices list.
        /// </summary>
        /// <value>The segment vertices list.</value>
        public List<Vertex> GetVertexes
        {
            get
            {
                return vertices;
            }
        }

        #endregion _properties

    }
}
