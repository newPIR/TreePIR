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
