#ifndef _FACE_H
#define _FACE_H


#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

struct Triangle {
    int node1, node2, node3;
    Triangle(int node1_in, int node2_in, int node3_in) : node1(node1_in), node2(node2_in), node3(node3_in) {}
    Triangle() : node1(0), node2(0), node3(0) {}
};

struct Node {
    double x, y, z;
    Node(double x_in, double y_in , double z_in) : x(x_in), y(y_in), z(z_in) {}
    Node() : x(0.0), y(0.0), z(0.0) {}
};

class Face {
    public:
        /**
         * Arg i;
         * xDimIn, yDimIn, zDimIn refers to thhe 'x/y/zDimNum', ie, number of grids on that dimension, the length of that dimension (x/y/zDim) should be x/y/zDimNum*delta
        */
        Face(int offsetIn, string fixedAxisIn, int fixedAxisValIn, int xDimIn, int yDimIn, int zDimIn, int deltaIn) : 
            offset(offsetIn), fixedAxis(fixedAxisIn), fixedAxisVal(fixedAxisValIn),
            xDimNum(xDimIn), yDimNum(yDimIn), zDimNum(zDimIn), delta(deltaIn),
            xDim(xDimNum* delta), yDim(yDimNum* delta), zDim(zDimNum* delta)
        {}
 
        inline size_t numNodes() { return nodes.size(); }
        inline size_t numTriangles() { return triangles.size(); }

        void genNodes() {
            int iLim, jLim;
            setLoopLimit(iLim, jLim, true);
            // Set iLim, jLim to be 'Dim' of two unfixed dimensions

            bool xIsFixed = fixedAxis.compare("x") == 0;
            bool yIsFixed = fixedAxis.compare("y") == 0;
            bool zIsFixed = fixedAxis.compare("z") == 0;

            for (int i = 0; i <= iLim; i += delta) {
                for (int j = 0; j <= jLim; j += delta) {
                    if (xIsFixed) {
                        nodes.emplace_back(fixedAxisVal, i, j + offset);
                    }
                    else if (yIsFixed) {
                        nodes.emplace_back(i, fixedAxisVal, j + offset);
                    }
                    else {
                        //? if z is fixed, meaning of offset?
                        nodes.emplace_back(i, j + offset, fixedAxisVal);
                    }
                }
            }
        }

        void shiftNodes(double deltaX, double deltaY, double deltaZ) {
            for (Node & n : nodes) {
                n.x += deltaX;
                n.y += deltaY;
                n.z += deltaZ;
            }
        }

        void genTriangles() {
            int iLim, jLim;
            setLoopLimit(iLim, jLim, false);

            for (int i = 0; i < iLim; ++i) {
                for (int j = 0; j < jLim; ++j) {
                    triangles.emplace_back(i * (jLim + 1) + j, (i * (jLim + 1)) + (j + 1), (i + 1) * (jLim + 1) + j);
                    triangles.emplace_back((i * (jLim + 1)) + (j + 1), (i + 1) * (jLim + 1) + j, (i + 1) * (jLim + 1) + (j + 1));
                }
            }
        }

        /**
         * @param isPoly, if yes, output in .poly form, else output in .off form
         * @param printTriangleHeader, is useful when isPoly is set true, if set to true, print a header to the triangle output..
        */
        void genOutput(int& nNodes, int& nTriangles, ostringstream& nodeOut, ostringstream& triangleOut, bool isPoly, bool printTriangleHeader) {
            if(isPoly){
                //.poly format output (triangle: index start from 1, node index printed..)
                if(printTriangleHeader){
                    triangleOut << "# A Facet..\n";
                    triangleOut << numTriangles() << endl;                    
                }
                for (int i = 0; i < triangles.size(); ++i) {
                    // nNodes HERE represent the offset (number of nodes already printed to the nodeOut string before this function call)
                    triangleOut << "3\t" << triangles[i].node1 + nNodes + 1
                        << "\t" << triangles[i].node2 + nNodes + 1
                        << "\t" << triangles[i].node3 + nNodes + 1 << "\n";
                    nTriangles++;
                }
                for (int i = 0; i < nodes.size(); ++i) {
                    nNodes++;
                    nodeOut << nNodes << "\t" << nodes[i].x << "\t" << nodes[i].y << "\t" << nodes[i].z << "\n";
                    }
            }else{
                //.off format output
                for (int i = 0; i < triangles.size(); ++i) {
                    // nNodes HERE represent the offset (number of nodes already printed to the nodeOut string before this function call)
                    triangleOut << "3\t" << triangles[i].node1 + nNodes
                        << "\t" << triangles[i].node2 + nNodes
                        << "\t" << triangles[i].node3 + nNodes << "\n";
                    nTriangles++;
                }
                for (int i = 0; i < nodes.size(); ++i) {
                    nodeOut << nodes[i].x << "\t" << nodes[i].y << "\t" << nodes[i].z << "\n";
                    nNodes++;
                    }
            }

        }

    private:
        vector<Node> nodes;
        vector<Triangle> triangles;
        string fixedAxis;
        int fixedAxisVal;
        int offset;
        int delta;
        int xDimNum, yDimNum, zDimNum;
        int xDim, yDim, zDim;

        void setLoopLimit(int& iLim, int& jLim, bool isNode) {
            if (isNode) {
                if (fixedAxis.compare("x") == 0) {
                    iLim = yDim;
                    jLim = zDim;
                }
                else if (fixedAxis.compare("y") == 0) {
                    iLim = xDim;
                    jLim = zDim;
                }
                else if (fixedAxis.compare("z") == 0) {
                    iLim = xDim;
                    jLim = yDim;
                }
                else {
                    cout << "Fixed Axis string malformed.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else {
                // Option == 0 => Triangles (xDimNum = xDim/delta)
                if (fixedAxis.compare("x") == 0) {
                    iLim = yDimNum;
                    jLim = zDimNum;
                }
                else if (fixedAxis.compare("y") == 0) {
                    iLim = xDimNum;
                    jLim = zDimNum;
                }
                else if (fixedAxis.compare("z") == 0) {
                    iLim = xDimNum;
                    jLim = yDimNum;
                }
                else {
                    cout << "Fixed Axis string malformed.\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
};




class NonuniformFace {
    public:

        // Note: by this definition, the delta given is for the middle part, while the DimNum is for the Upper and Bottom part.
        NonuniformFace(string fixedAxisIn, int fixedAxisValIn, int xDimNumIn, int yDimNumIn, int zDimNumTIn, int zDimNumBIn, int zDimNumMIn, int offsetIn, int deltaIn) : 
            bottom(0, fixedAxisIn, fixedAxisValIn, xDimNumIn, yDimNumIn, zDimNumBIn, 2 * deltaIn),
            middle(offsetIn + zDimNumBIn * 2 * deltaIn, fixedAxisIn, fixedAxisValIn, xDimNumIn * 2, yDimNumIn * 2, zDimNumMIn, deltaIn),
            top(2 * offsetIn + zDimNumMIn * deltaIn + zDimNumBIn * 2 * deltaIn, fixedAxisIn, fixedAxisValIn, xDimNumIn, yDimNumIn, zDimNumTIn, 2 * deltaIn),
            fixed_axis(fixedAxisIn), xDimNum(xDimNumIn), yDimNum(yDimNumIn), zDimTNum(zDimNumTIn),zDimBNum(zDimNumBIn), zDimMNum(zDimNumMIn), delta(deltaIn)
        {}

        void genNodes() {
            bottom.genNodes();
            middle.genNodes();
            top.genNodes();
        }

        void shiftNodes(double deltaX, double deltaY, double deltaZ) {
            bottom.shiftNodes(deltaX, deltaY, deltaZ);
            middle.shiftNodes(deltaX, deltaY, deltaZ);
            top.shiftNodes(deltaX, deltaY, deltaZ);
        }

        void genTriangles() {
            bottom.genTriangles();
            middle.genTriangles();
            top.genTriangles();
        }

        void genConnectorsBottom(int dim) {
            int startingIndexBottom = zDimBNum;
            int startingIndexMiddle = bottom.numNodes();
            for (int i = 0; i < dim; i++) {
                connectingLayer.emplace_back(startingIndexBottom + i * (zDimBNum + 1), startingIndexMiddle + 2 * i * (zDimMNum + 1), startingIndexMiddle + (2 * i + 1) * (zDimMNum + 1));
                connectingLayer.emplace_back(startingIndexBottom + i * (zDimBNum + 1), startingIndexBottom + (i + 1) * (zDimBNum + 1), startingIndexMiddle + (2 * i + 1) * (zDimMNum + 1));
                connectingLayer.emplace_back(startingIndexBottom + (i + 1) * (zDimBNum + 1), startingIndexMiddle + (2 * i + 1) * (zDimMNum + 1), startingIndexMiddle + 2 * (i + 1) * (zDimMNum + 1));
            }
        }

        void genConnectorsTop(int dim) {
            size_t numNodesBottomFace = bottom.numNodes();
            size_t numNodesMiddleFace = middle.numNodes();
            int startingIndexTop = numNodesBottomFace + numNodesMiddleFace;
            int startingIndexMiddle = numNodesBottomFace + zDimMNum;
            for (int i = 0; i < dim; i++) {
                connectingLayer.emplace_back(startingIndexTop + i * (zDimTNum + 1), startingIndexMiddle + 2 * i * (zDimMNum + 1), startingIndexMiddle + (2 * i + 1) * (zDimMNum + 1));
                connectingLayer.emplace_back(startingIndexTop + i * (zDimTNum + 1), startingIndexTop + (i + 1) * (zDimTNum + 1), startingIndexMiddle + (2 * i + 1) * (zDimMNum + 1));
                connectingLayer.emplace_back(startingIndexTop + (i + 1) * (zDimTNum + 1), startingIndexMiddle + (2 * i + 1) * (zDimMNum + 1), startingIndexMiddle + 2 * (i + 1) * (zDimMNum + 1));
            }
        }


        void genTrianglesAndConnectors(int dimNumWidth){
            genTriangles();
            genConnectorsBottom(dimNumWidth);
            genConnectorsTop(dimNumWidth);
        }

        /**
         * @param isPoly, if yes, output in .poly form, else output in .off form
        */
        void printConnectors(int & nNodes, int & nTriangles, ostringstream & triangleOut, vector<Triangle> & toPrint, bool isTop, bool isPoly){
            int nodeCorrection = bottom.numNodes();
            if (isTop) {
                nodeCorrection += (middle.numNodes() + top.numNodes());
            }
            if (isPoly) {
                nodeCorrection -= 1;
            }
            for (int i = 0; i < toPrint.size(); ++i) {
                triangleOut << "3\t" << toPrint[i].node1 + nNodes - nodeCorrection
                            << "\t" << toPrint[i].node2 + nNodes - nodeCorrection
                            << "\t"  << toPrint[i].node3 + nNodes - nodeCorrection << "\n";
                nTriangles++; 
            } 
        }

        /**
         * @param isPoly, if yes, output in .poly form, else output in .off form
         * @note in .off mode, this won't print the number of triangles in each facet into the output, need to do this elsewhere
        */
        void genOutput(int & nNodes, int & nTriangles, ostringstream & nodeOut, ostringstream & triangleOut, bool isPoly) {

            triangleOut << "# A Non-uniform Facet..\n";
            triangleOut << numTriangles() << endl;
            // triangleOut << "BOTTOM Triangles START HERE\n";
            bottom.genOutput(nNodes, nTriangles, nodeOut, triangleOut, isPoly, false);
            // triangleOut << "BOTTOM Triangles END HERE\n";
            //triangleOut << "### First Connector" << endl;
            printConnectors(nNodes, nTriangles, triangleOut, connectingLayer, false, isPoly);
            //triangleOut << "### END: First Connector" << endl;
            // triangleOut << "MIDDLE Triangles START HERE\n";
            middle.genOutput(nNodes, nTriangles, nodeOut, triangleOut, isPoly, false);
            // triangleOut << "TOP Triangles START HERE\n";
            top.genOutput(nNodes, nTriangles, nodeOut, triangleOut, isPoly, false);
            //triangleOut << "### Second Connector" << endl;
            //printConnectors(nNodes, nTriangles, triangleOut, connectingLayer, true, isPoly);

        }

        inline int numTriangles() {
            return top.numTriangles() + middle.numTriangles() + bottom.numTriangles() + connectingLayer.size();
        }

    private:
        Face top;
        Face middle;
        Face bottom;
        string fixed_axis;
        int xDimNum, yDimNum, zDimTNum, zDimBNum, zDimMNum, delta;
        vector<Triangle> connectingLayer;
};

#endif