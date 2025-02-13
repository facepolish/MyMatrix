// swift-tools-version: 6.0
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "MyMatrix",
    platforms: [
        .macOS(.v12), // 例: macOS 12 以降
        .iOS(.v13),   // 例: iOS 13 以降
        // 必要に応じて他のプラットフォームも追加
    ],
    products: [
        .library(name: "MyMatrix", targets: ["MyMatrix"]), // ライブラリ製品を定義
        .executable(name: "MyMatrixExe", targets: ["MyMatrixExe"]), // 実行可能ファイル製品を定義 (必要であれば)
    ],
    targets: [
        .target(name: "MyMatrix"), // ライブラリターゲット
        .executableTarget(name: "MyMatrixExe", dependencies: ["MyMatrix"]), // 実行可能ファイルターゲット (必要であれば)
    ]
)
