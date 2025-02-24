// The Swift Programming Language
// https://docs.swift.org/swift-book
import Foundation
import MyMatrix
/*
func sukkiri(_ a:Matrix)->Matrix {
    let epsilon:Float = 0.00001
    var b:[Float] = []
    b = a.flat.map{ check($0)}
    return Matrix(b,row:a.row, col:a.col)
    func check(_ val:Float) -> Float {
        if (val < epsilon) {
            return 0
        }else{
            return val
        }
    }
}

let a1: [[Float]] = [
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [9.0, 10.0, 11.0, 12.0],
    [13.0, 14.0, 15.0, 16.0]
]

let b1: [[Float]] = [
    [17.0, 18.0, 19.0, 20.0],
    [21.0, 22.0, 23.0, 24.0],
    [25.0, 26.0, 27.0, 28.0],
    [29.0, 30.0, 31.0, 32.0]
]
let sample:[[Float]] = [
    [1,1,1,-1],
    [1,1,-1,1],
    [1,-1,1,1],
    [-1,1,1,1]
]
let c1:[[Float]] = [
    [1,3,0,-2],
    [4,0,3,5],
    [5,4,2,0],
    [0,1,3,2]
]
do {
    let a2 = try randMatrix(4).orthnormal()
    let b2 = a2.transpose()
    let c2 = try b2 * a2
    print(sukkiri(c2).get())
    
}
func makeClear(_ a:Matrix)->Matrix {
    let epsilon:Float = 0.00001
    var b:[Float] = []
    b = a.flat.map{ check($0)}
    return Matrix(b,row:a.row, col:a.col)
    func check(_ val:Float) -> Float {
        if (val < epsilon) {
            return 0
        }else{
            return val
        }
    }
}
func compare(_ a: [[Float]],_ b:[[Float]]) -> Bool {
    let col = a.count
    let row = a[0].count
    guard a.count == b.count && a[0].count == b[0].count else {
        return false
    }
    for i in 0..<row {
        for j in 0..<col {
            if abs(a[i][j] - b[i][j]) > 0.01 {
                return false
            }
        }
    }
    return true
}

let checkMatrix1:[[Float]] = [[1,0,0],[0,1,0],[0,0,1]]
let checkMatrix2:[[Float]] = [[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]]
do {
    let matrix1 = try randMatrix(3)
    print(matrix1)
    let a1 = try matrix1.orthnormal()
    let b1 = a1.transpose()
    let c1 = try b1 * a1
    print(makeClear(c1).get())
    if compare(c1.get(),checkMatrix1){
        print ("test is OK")
    }else{
        print("test is NG")
    }
    
    let matrix2 = try randMatrix(5)
    print(matrix2)
    let a2 = try matrix2.orthnormal()
    let b2 = a2.transpose()
    let c2 = try b2 * a2
    print(makeClear(c2).get())
    if compare(c2.get(),checkMatrix2){
        print ("test is OK")
    }else{
        print("test is NG")
    }
}
*/
let epsilon:Float = 0.0001
func printEpsilon() {
    print (epsilon)
}
