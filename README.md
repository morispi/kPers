# kPers

## Build

```
mkdir build
cmake ..
make -j4
```

If cmake didn't find your seqan install you can run

```
mkdir build
cmake -DCMAKE_CXX_FLAGS="-I /path/to/your/seqan/include" ..
make -j4
```
