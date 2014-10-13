//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#include "diophantine.h"

const std::map<dioph::id,iArray> dioph::nodes =
{
    // dim = 2
    { dioph::d2n52,   {2,   52,   15,   25} },
    { dioph::d2n538,  {2,  538,  171,  177} },
    { dioph::d2n1154, {2, 1154,  177,  415} },
    { dioph::d2n3722, {2, 3722, 1259, 1403} },
    { dioph::d2n6044, {2, 6044, 1427, 1891} },
    
    // dim = 3
    { dioph::d3n52,   {3,   52,    5,    7,   11} },
    { dioph::d3n538,  {3,  538,  115,  241,  251} },
    { dioph::d3n1154, {3, 1154,   11,  369,  569} },
    { dioph::d3n3722, {3, 3722,  653, 1005, 1483} },
    { dioph::d3n6044, {3, 6044,  595, 1359, 1865} },
    
    // dim = 4
    { dioph::d4n1154, {4, 1154, 103,  217,  453,  521} },
    { dioph::d4n3722, {4, 3722, 119,  203,  235, 1849} },
    { dioph::d4n6044, {4, 6044, 687, 1385, 2843, 2933} },
    
    // dim = 5
    { dioph::d5n1154, {5, 1154, 107,  211,  225,  459,  547} },
    { dioph::d5n3722, {5, 3722, 263, 1049, 1145, 1441, 1759} },
    { dioph::d5n6044, {5, 6044, 681, 1113, 1165, 1853, 2349} },
    
    // dim = 6
    { dioph::d6n1154, {6, 1154,   17,   71,  151,  203,  209,  385} },
    { dioph::d6n2008, {6, 2008,  177,  451,  453,  565,  665,  833} },
    { dioph::d6n3722, {6, 3722,  205,  399,  679, 1051, 1163, 1811} },
    { dioph::d6n6044, {6, 6044,  587,  671,  673, 1459, 2363, 2697} },
    { dioph::d6n9644, {6, 9644, 1177, 1269, 1625, 1947, 3369, 4087} },
    
    // dim = 7
    { dioph::d7n1154, {7, 1154,   3,   5,   37,  377,  393,  479,  559} },
    { dioph::d7n3722, {7, 3722, 669, 819,  845,  857, 1379, 1421, 1439} },
    { dioph::d7n6044, {7, 6044,  35, 397, 1427, 1467, 1891, 2805, 2865} },
    
    // dim = 8
    { dioph::d8n3708, {8, 3708, 151, 187,  333,  357,  839,  947, 1245, 1839} },
    { dioph::d8n6044, {8, 6044,  43,  61,  179, 1421, 1479, 2039, 2189, 2783} },
    { dioph::d8n9644, {8, 9644, 129, 525, 2373, 2731, 3351, 4145, 4419, 4767} },
    
    // dim = 9
    { dioph::d9n3722, {8, 3722, 119, 339,  437,  773,  937, 1219, 1503, 1697, 1747} },
    { dioph::d9n6044, {8, 6044,  43,  87,  179, 1421, 1479, 1589, 2189, 2191, 2783} },
    { dioph::d9n9644, {8, 9644, 457, 509, 1677, 1723, 2173, 2423, 2489, 3431, 3719} },
    
    // dim = 10
    { dioph::d10n3722, {10, 3722, 153, 223, 517,  859,  861,  911,  991,  995, 1453, 1705} },
    { dioph::d10n6044, {10, 6044,  43,  87,  89,  179, 1421, 1479, 1589, 2095, 2191, 2783} },
    { dioph::d10n9644, {10, 9644,  53, 263, 893, 1121, 1259, 2945, 3281, 4397, 4727, 4801} },
    
    // dim = 11
    { dioph::d11n3722, {11, 3722,  39,  41,  51,  481,  649,  671,  705, 1213, 1297, 1559, 1719} },
    { dioph::d11n6044, {11, 6044, 159, 289, 309,  349,  471,  523,  727,  847,  951, 1665, 2899} },
    { dioph::d11n9644, {11, 9644, 205, 723, 731, 1693, 2215, 2383, 2583, 2787, 2961, 3061, 3123} },
    
    // dim = 12
    { dioph::d12n3722, {12, 3722,  97,  419,  435,  477,  641,  717,  793, 1141, 1163, 1359, 1535, 1787} },
    { dioph::d12n6044, {12, 6044, 179,  795, 1073, 1233, 1359, 1589, 1667, 1675, 1677, 1917, 2095, 2189} },
    { dioph::d12n9644, {12, 9644, 973, 1095, 1495, 1853, 2619, 3365, 3493, 3567, 3595, 3863, 4063, 4497} }
};
