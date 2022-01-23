from typing import Optional, List, Tuple, Any

import matplotlib.pyplot as plt

from .scoring_matrix import PKModule


def draw_logo(
    state_probs: List[List[Tuple[Any, float]]],
    save_to: Optional[str] = None
) -> None:

    x = [i for i in range(len(state_probs))]
    names = [i + 1 for i in x]

    for module_idx, module in enumerate(PKModule):
        scores = [
            state_prob[module_idx][1] for state_prob in state_probs
        ]
        plt.bar(
            x,
            scores,
            color=module.logo_color(),
            edgecolor='white',
            width=1.0
        )

    plt.xticks(x, names, fontweight='bold')
    plt.xlabel('position')
    plt.ylabel('probability')
    plt.title('conserved profile')

    labels = [module.display_name() for module in PKModule][:-1]  # Exclude gap
    colors = [module.logo_color() for module in PKModule][:-1]  # Exclude gap
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors[idx])
        for idx, _ in enumerate(labels)
    ]
    plt.legend(handles, labels, bbox_to_anchor=(1.2, 0.925))
    plt.tight_layout()

    if save_to is None:
        plt.show()
        plt.clf()
    else:
        plt.savefig(save_to, dpi=300)


