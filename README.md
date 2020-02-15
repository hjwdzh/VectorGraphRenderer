# VectorGraphRenderer
A renderer that takes a triangle mesh, a camera pose and produce a vector graph.

![VectorGraphRender Teaser](https://github.com/hjwdzh/VectorGraphRenderer/raw/master/res/teaser.png)

### Compile
```
mkdir build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```

### Run example
Run binary with three arguments. The first is the input mesh, the second the camera parameters, and the third the output file.

If the output file ends with svg, we produce a svg image file. If the output file ends with obj, we produce a obj mesh output.

Example of running the script:
```
./vector_render ../example/AccuCities_HD.obj ../example/matrix.txt result.svg
```

