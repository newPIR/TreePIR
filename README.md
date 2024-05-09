<p align="center">
  <img width="250" height="230" src="https://github.com/PIR-PIXR/Certificate-Transparency-Logs/assets/102839948/530caacf-868e-464c-995b-04e995bc02bc">
</p>

# TreePIR: Efficient Private Retrieval of Merkle Proofs via Tree Colorings with Fast Indexing and Zero Storage Overhead

## Abstract
A Batch Private Information Retrieval (batch-PIR) scheme allows a client to retrieve multiple data items from a database without revealing them to the storage server(s). Most existing approaches for batch-PIR are based on batch codes, which incur large storage overheads. In this work, we show that *zero* storage overhead is achievable for tree-shaped databases. In particular, we develop *TreePIR*, a novel approach tailored made for private retrieval of the set of nodes along an arbitrary *root-to-leaf path* in a Merkle tree with no storage redundancy. This type of trees has been widely implemented in many real-world systems such as Amazon DynamoDB, Google's Certificate Transparency, and blockchains. Tree nodes along a root-to-leaf path forms the well-known *Merkle proof*. TreePIR, which employs a novel tree coloring, outperforms the state-of-the-art Probabilistic Batch Codes (PBC) ([Angel et al. IEEE S&P'18](https://eprint.iacr.org/2017/1142.pdf)) in all metrics, achieving $3\times$ lower total storage and communication cost and $1.5$-$2\times$ faster max/total server computation time and client query generation time. Most notably, TreePIR's *polylog*-complexity indexing algorithm is $19$ - $160\times$ faster than PBC for trees of $2^{10}$ - $2^{24}$ leaves.

---
## Experimental setup
Using only one core, our experiments ran within Ubuntu 22.04 LTS environments (Intel® Core™ i9-13900H and 32GiB of system memory). Each experiment was repeated ten times for various trees. The average values were then calculated. We fetched $2^{20}$ entries from Google's [Xenon2024](https://github.com/PIR-PIXR/Certificate-Transparency-Logs) for running end-to-end PIR system. Each entry comprised an entry number, timestamp, and certificate. Applying SHA-256 on the certificates, we constructed Merkle trees of $n$ leaves with $n$ ranging from $2^{10}$ to $2^{20}$. Consequently, each tree node occupied 32 bytes.

To benchmark TreePIR and PBC for large trees ($n = 2^{22},\ldots,2^{36})$, we use random hashes to avoid excessive overheads.


---
## Installing Libraries

- #### Javac
      $ sudo apt update
      $ sudo apt upgrade
      $ sudo apt install default-jdk
- #### SEAL 4.0.0
      $ sudo apt install build-essential cmake clang git g++ libssl-dev libgmp3-dev
      $ sudo apt update
      $ sudo apt upgrade
      $ git clone https://github.com/cnquang/SEAL-4.0.0.git
      $ cd SEAL-4.0.0
      $ cmake -S . -B build
      $ cmake --build build
      $ sudo cmake --install build
- #### JSON
      $ git clone https://github.com/microsoft/vcpkg
      $ ./vcpkg/bootstrap-vcpkg.sh
      $ ./vcpkg install rapidjson
- #### Google gRPC
      $ sudo apt install -y build-essential autoconf libtool pkg-config
      $ git clone --recurse-submodules -b v1.58.0 --depth 1 --shallow-submodules https://github.com/grpc/grpc
      $ cd grpc
      $ mkdir -p cmake/build
      $ pushd cmake/build
      $ cmake -DgRPC_INSTALL=ON \
        -DgRPC_BUILD_TESTS=OFF \
        ../..
      $ make -j 4
      $ sudo make install
      $ popd

---
### Executing SealPIR+PBC and SealPIR+TreePIR
      $ git clone https://github.com/newPIR/TreePIR.git
      $ cd TreePIR/SealPIR-Orchestrator
      $ python3 orchestrator.py <tree height> /Path/To/TreePIR
For example, with tree height $h = 10$, 
      python3 orchestrator.py 10 /Path/To/TreePIR.

---
### Executing Spiral+PBC and Spiral+TreePIR

- ##### Docker
      $ cd Spiral
      $ sudo docker build -t spiral_toolchain .
      $ sudo apt-get install libstdc++-11-dev
      $ sudo docker run -it \
        -u root \
        -v /Path/To/TreePIR/Spiral:/tmp/Spiral \
        -v /Path/To/TreePIR/Spiral/Process_Workspace:/home/ubuntu/Process_Workspace \
        --rm spiral_toolchain:latest \
        /bin/bash -c "cd /tmp/Spiral; exec bash"
      $ cd Evaluation
      $ python3 evaluate.py

  For executing Spiral+PBC, you should open the evaluation.py and active line 23: "mode: typing.Final[Mode] = Mode.PBC". And in the file TreePIR/Spiral/Seperate/src/spiral.cpp, active line 1140.

  For executing Spiral+TreePIR, you should open the evaluation.py and active line 22: "mode: typing.Final[Mode] = Mode.COLOURING". And in the file TreePIR/Spiral/Seperate/src/spiral.cpp, active line 1139.





      
