from __future__ import annotations
from enum import Enum, auto


import matplotlib.pyplot as plt


class PKModule(Enum):
    A1, A2, A3 = auto(), auto(), auto()
    B1, B2 = auto(), auto()
    C = auto()
    D = auto()
    E1, E2 = auto(), auto()
    F1, F2 = auto(), auto()
    G1, G2 = auto(), auto()
    H = auto()
    I1, I2 = auto(), auto()
    GAP = auto()  # Denotes gap in PK module sequence

    def display_name(self):
        if self == PKModule.GAP:
            return '-'
        else:
            return self.name

    def display_in_alignment(self):
        max_length = max([len(m.display_name()) for m in PKModule])
        padding = max_length - len(self.display_name())
        return (' ' * padding) + self.display_name()

    def logo_color(self):
        cm = plt.get_cmap('gist_rainbow')
        colors = {
            module.display_name(): cm(1. * module_idx/len(PKModule))
            for module_idx, module in enumerate(PKModule)
        }
        return colors[self.display_name()]


class ScoringMatrix:
    # Default scoring matrix:
    #   score  2 if modules are identical
    #   score -1 if they are non-identical
    #   score -4 either module is a gap

    # Build default scoring matrix:
    pairwise_scores = {
        module1: {
            module2: (
                2 if module1 == module2
                else -4 if (
                    module1 == PKModule.GAP or
                    module2 == PKModule.GAP
                )
                else -1
            )
            for module2 in PKModule
        }
        for module1 in PKModule
    }

    def __call__(self, module1: PKModule, module2: PKModule) -> float:
        try:  # Validate scoring matrix on first call
            if self.__class__.__call__._called:
                return self._score(module1, module2)  # Score modules
        except AttributeError:
            if not self._validate():  # Validate scoring matrix
                raise ValueError('ScoringMatrix is invalid')
            self.__class__.__call__._called = True
            # Call scoring matrix again after validation
            return self(module1, module2)

    def _validate(self) -> bool:
        # Validate if all modules are pairwise represented in scoring matrix
         return all([
            all([
                m2 in self.pairwise_scores[m1]
                for m2 in PKModule
            ])
            for m1 in PKModule
        ])

    def _score(self, module1: PKModule, module2: PKModule) -> float:
        return self.pairwise_scores[module1][module2]
