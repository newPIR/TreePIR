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

To execute Spiral+PBC, open the evaluation.py and active line 23: "mode: typing.Final[Mode] = Mode.PBC". In the file TreePIR/Spiral/Seperate/src/spiral.cpp, active line 1140,

To execute Spiral+TreePIR, open the evaluation.py and active line 22: "mode: typing.Final[Mode] = Mode.COLOURING". In the file TreePIR/Spiral/Seperate/src/spiral.cpp, active line 1139,

---
### Executing VBPIR+PBC
      $ cd VBPIR_PBC
      $ cmake -S . -B build
      $ cmake --build build
      $ ./build/bin/vbpir_pbc
      
---
### Executing VBPIR+TreePIR
      $ cd VBPIR_TreePIR
      $ cmake -S . -B build
      $ cmake --build build
      $ ./build/bin/vbpir_treepir

---
## Performance
The performance of a batch-code-based batch-PIR depends on the number of sub-databases and their sizes, and on the performance of the underlying PIR scheme. Moreover, it also depends on the dimension $d$ that PIR sub-databases are represented, e.g., $d=2$ for SealPIR, $d=3$ for VBPIR, and $d=4$ for Spiral. Theoretically, as TreePIR uses $1.5\times$ fewer sub-databases of size $2\times$ smaller compared to PBC, theoretically, TreePIR uses $3\times$ less storage, with $\sqrt[d]{2}\times$ faster max server computation and client query-generation times, $1.5\sqrt[d]{2}\times$ faster total server computation time. The client answer-extraction times are similar for both as dummy responses are ignored in PBC. The theoretical gains were also reflected in the experiments.

TreePIR has significantly faster setup and indexing thanks to its efficient coloring and indexing algorithms on trees. In particular, for trees with $2^{10}$ - $2^{24}$ leaves, TreePIR's setup and indexing are $8$-$60\times$ and $19$-$160\times$ faster than PBC's, respectively. TreePIR still works well beyond that range, requiring $180$ seconds to setup a tree of $2^{30}$ leaves, and only $0.7$ milliseconds to index in a tree of $2^{36}$ leaves.

| $h$ | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 |
|-----|----|----|----|----|----|----|----|----|
| PBC (ms) | 3.4 | 8.9 | 53.1 | 406 | 2132 | 9259 | 37591 | 159814 |
| **TreePIR** (ms) | 0.4 | 2.4 | 7.6 | 39.1 | 56.1 | 179 | 600 | 2700 |
| $h$ | 26 |    |    | 28 |    | 29 |    | 30 |
| **TreePIR** (sec) | 9.6 |    |    | 37.3 |    | 77.1 |    | 179.4 |

Table 1: A comparison of the *setup* time between TreePIR and PBC for various tree heights. TreePIR's setup is *8*-*$60\times$ faster.

<p align="center">
  <img width="400" height="300" src="https://user-images.githubusercontent.com/87842051/230779728-d474a1e3-6d0a-4ad2-83d4-5ab550d265fa.png"> <img width="400" height="300" src="https://user-images.githubusercontent.com/87842051/230779764-8e25f985-b584-4f24-9cab-def75f1c7efa.png">
</p>
<strong> Fig. 2.</strong> The comparison of the average server and client computation times between LM-CKGS and LM-BE schemes where n = 2^10, m = 2^12, t = (k − 1) for both schemes, t = 1 for LM-BE scheme, and k increases from 3 to 6. Assume (k − 1) colluding servers in the system, the LM-BE scheme performs quite better on the server side. In contrast, the LM-CKGS computational cost is slightly lower on the Client side. For the LM-BE scheme, there is a tradeoff when t changes and higher t means higher computation costs.

<p align="center">
  <img width="400" height="300" src="https://user-images.githubusercontent.com/87842051/230779833-3170bd8e-fe30-4c38-9a2c-93dcd0349aad.png"> <img width="400" height="300" src="https://user-images.githubusercontent.com/87842051/230779844-739ba227-cd21-4eb8-8732-4354c653d709.png">
</p>
<strong> Fig. 3.</strong> The comparison of the average server and client computation times regraded to three different LM-PIR cases (LM-CKGS, LM-WY, and LM-BE) where k = 2, t = 1, m = 2^10, and n increases from 2^8 to 2^12. When n increases, the total LM-PIR computation time increases, especially the LM computation time increases significantly. It means that the percentage of LMC-related computation time increases significantly when the horizontal database increases.

<p align="center">
  <img width="400" height="300" src="https://user-images.githubusercontent.com/87842051/230779906-eb2ffdf5-a405-4989-9ee5-f0d9c51530b7.png"> <img width="400" height="300" src="https://user-images.githubusercontent.com/87842051/230779912-8aeb452f-e90f-4947-978c-4230d8fddf64.png">
</p>
<strong> Fig. 4.</strong> The comparison of the average server and client computation times of the LM-WY scheme where t = 1, n = 2^10, m = (2^10, 2^12), and (d, k) =
((3, 2), (4, 3), (5, 3), (6, 4)). When d increases, the LM-PIR computation time on the server side tends to decrease because the number of answers and witnesses decreases regarded to O(n^(1/d)). However, when d increases, the size of queries also grows. It is why the computation cost of higher d is slightly higher in some cases, but in general, the computation time on the server side reduces when d rises. On the Client side, the computation cost trend is similar to the computation time on the server. However, d increases lead to an increase in k, so the total computation time on the Client side grows.

---
## REFERENCES

[LM] R. W. Lai and G. Malavolta, “Subvector commitments with application to succinct arguments,” in Advances in Cryptology–CRYPTO 2019: 39th Annual International Cryptology Conference, Santa Barbara, CA, USA, August 18–22, 2019, Proceedings, Part I 39. Springer, 2019, pp. 530–560.

[CKGS] B. Chor, E. Kushilevitz, O. Goldreich, and M. Sudan, “Private information retrieval,” Journal of the ACM (JACM), vol. 45, no. 6, pp. 965–981, 1998.

[WY] D. Woodruff and S. Yekhanin, “A geometric approach to information-theoretic private information retrieval,” in 20th Annual IEEE Conference on Computational Complexity (CCC’05). IEEE, 2005, pp. 275–284.

[BE] R. Bitar and S. El Rouayheb, “Staircase-pir: Universally robust private information retrieval,” in 2018 IEEE Information Theory Workshop (ITW). IEEE, 2018, pp. 1–5.
