from typing import Optional, Sequence, Tuple, Union

import cudf
import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from packaging import version
from scanpy.tools._utils import _choose_graph
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency
from scipy import sparse


def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    key_added: str = "leiden",
    adjacency: Optional[sparse.spmatrix] = None,
    n_iterations: int = 100,
    use_weights: bool = True,
    neighbors_key: Optional[str] = None,
    obsp: Optional[str] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Performs Leiden clustering using cuGraph, which implements the method
    described in:

    Traag, V.A., Waltman, L., & van Eck, N.J. (2019). From Louvain to
    Leiden: guaranteeing well-connected communities. Sci. Rep., 9(1), 5233.
    DOI: 10.1038/s41598-019-41695-z

    Parameters
    ----------
        adata :
            annData object

        resolution
            A parameter value controlling the coarseness of the clustering.
            (called gamma in the modularity formula). Higher values lead to
            more clusters.

        random_state
            Change the initialization of the optimization. Defaults to 0.

        restrict_to
            Restrict the clustering to the categories within the key for
            sample annotation, tuple needs to contain
            `(obs_key, list_of_categories)`.

        key_added
            `adata.obs` key under which to add the cluster labels.

        adjacency
            Sparse adjacency matrix of the graph, defaults to neighbors
            connectivities.

        n_iterations
            This controls the maximum number of levels/iterations of the
            Leiden algorithm. When specified, the algorithm will terminate
            after no more than the specified number of iterations. No error
            occurs when the algorithm terminates early in this manner.

        use_weights
            If `True`, edge weights from the graph are used in the
            computation (placing more emphasis on stronger edges).

        neighbors_key
            If not specified, `leiden` looks at `.obsp['connectivities']`
            for neighbors connectivities. If specified, `leiden` looks at
            `.obsp[.uns[neighbors_key]['connectivities_key']]` for neighbors
            connectivities.

        obsp
            Use .obsp[obsp] as adjacency. You can't specify both
            `obsp` and `neighbors_key` at the same time.

        copy
            Whether to copy `adata` or modify it in place.
    """
    # Adjacency graph
    from cugraph import Graph
    from cugraph import leiden as culeiden

    adata = adata.copy() if copy else adata

    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories,
            adjacency,
        )
    offsets = cudf.Series(adjacency.indptr)
    indices = cudf.Series(adjacency.indices)
    if use_weights:
        weights = cudf.Series(adjacency.data)
    else:
        weights = None

    g = Graph()

    g.from_cudf_adjlist(offsets, indices, weights)

    # Cluster
    leiden_parts, _ = culeiden(
        g,
        resolution=resolution,
        max_iter=n_iterations,
    )

    # Format output
    groups = (
        leiden_parts.to_pandas().sort_values("vertex")[["partition"]].to_numpy().ravel()
    )
    if restrict_to is not None:
        if key_added == "leiden":
            key_added += "_R"
        groups = rename_groups(
            adata,
            key_added,
            restrict_key,
            restrict_categories,
            restrict_indices,
            groups,
        )
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype("U"),
        categories=natsorted(map(str, np.unique(groups))),
    )
    # store information on the clustering parameters
    adata.uns["leiden"] = {}
    adata.uns["leiden"]["params"] = {
        "resolution": resolution,
        "n_iterations": n_iterations,
    }
    return adata if copy else None