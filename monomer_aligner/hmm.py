from __future__ import annotations
from typing import Optional, Union, Dict, List
from enum import Enum, auto
from random import random
from math import log

from .parser import parse_fasta
from .scoring_matrix import PKModule
from .aligner import ModuleSequence, Module
from .drawer import draw_logo


# Background polyketide module probabilities:
background_probs = {
    module: 1 / (len(PKModule) - 1)  # Minus one for the GAP
    for module in PKModule
}


def normalize(counts: Dict[str, int]) -> Dict[str, float]:
    tot = sum([count for (key, count) in counts.items()])
    if tot:
        for key in counts:
            counts[key] /= tot
    return counts


def sample(
    events: Dict[Union[str, int, float], float]
) -> Union[str, int, float]:
    r = random()
    cum = 0
    for key, prob in events.items():
        cum += prob
        if r <= cum:
            return key


def logit(prob: float) -> float:
    try:
        calculated_logit = log(prob / (1 - prob))
    except ZeroDivisionError:
        raise ValueError('did you try to take log-odds of zero? Make sure'
                         'you use pseudo counts to prevent zeros')
    return calculated_logit


class HMM:
    def __init__(
        self,
        match_states: List[bool],
        pseudo_count: Optional[float] = 0.01
    ) -> None:
        self.match_states = match_states
        self.nmatches = sum(match_states)
        init_count = pseudo_count

        # Emission probabilities for match and insert states:
        self.e_m = {
            n: {
                a: init_count
                for a in background_probs
            }
            for n in range(self.nmatches)
        }
        self.e_i = background_probs

        # Transition probabilities from/to matches, inserts and deletions:
        self.t = {
            nmatch: {
                state: {
                    state: init_count
                    for state in ['m', 'i', 'd']
                }
                for state in ['m', 'i', 'd']
            }
            for nmatch in range(self.nmatches + 1)
        }

    def train(self, msa: List[ModuleSequence]) -> None:
        for seq in msa:
            self.count_states(seq)
        self.normalize_states()

    def count_states(self, seq: ModuleSequence) -> None:
        nmatch = 0
        route = [State(HMMStartEnd.START)]  # Sets start state

        for idx, ((module, _), match) \
        in enumerate(zip(seq._seq, self.match_states)):
            # Increments match state count after visiting match state:
            if match:
                nmatch = route[idx].nmatch + 1
                # Store emission state count:
                if module != PKModule.GAP:
                    self.e_m[nmatch - 1][module] += 1

            # Determine next state (m/i/d) based on previous state:
            state = State(module, match, nmatch)
            state.determine_state(route[idx])
            route.append(state)

        # Sets end state:
        route.append(State(HMMStartEnd.END, nmatch=(self.nmatches + 1)))

        # Filter states in route:
        route = list(
            filter(
                lambda state: state.module != PKModule.GAP
                or (not state.module == PKModule.GAP and state.match),
                route
            )
        )

        # Store transition state counts from filtered route:
        for s in range(len(route) - 1):
            b_s = route[s].state
            f_s = route[s + 1].state
            idx = route[s].nmatch
            self.t[idx][b_s][f_s] += 1

    def normalize_states(self) -> None:
        for nmatch in range(self.nmatches + 1):
            for start in self.t[nmatch]:
                self.t[nmatch][start] = normalize(self.t[nmatch][start])
            if nmatch < self.nmatches:
                self.e_m[nmatch] = normalize(self.e_m[nmatch])

    def sample_sequence(self) -> List[Module]:
        tag = None
        nmatch = 0
        emissions = []
        current = 'm'  # Begin state
        while nmatch < self.nmatches:
            next = sample(self.t[nmatch][current])

            # It can happen that training sequences do not contain certain
            # transitions (e.g. when number of sequences are low sequences
            # are not diverse enough). In that case, restart function:
            if next is None:
                try:
                    return self.sample_sequence()
                except RecursionError as err:
                    print(f'{err}: HMM contains too many unseen states in able'
                          f'to robustly sample sequences')

            current = next
            if current == 'm':
                emissions.append((sample(self.e_m[nmatch]), tag))
                nmatch += 1
            elif current == 'd':
                nmatch += 1
            elif current == 'i':
                emissions.append((sample(self.e_i), tag))
        return emissions

    def logo(self, save_to: Optional[str] = None) -> None:
        nmatch = 0
        state_logits = []
        for idx, match_state in enumerate(self.match_states):
            if match_state:
                probs = []
                for module in PKModule:
                    # prob = logit(self.e_m[nmatch][module])
                    # prob = 0.0 if prob < 0.0 else prob
                    prob = self.e_m[nmatch][module]
                    probs.append((module, prob))
                state_logits.append(probs)
                nmatch += 1
        draw_logo(state_logits, save_to)


class HMMStartEnd(Enum):
    START = auto()
    END = auto()


class State:
    def __init__(
        self,
        module: Union[PKModule, HMMStartEnd],
        match: Optional[bool] = False,
        nmatch: Optional[int] = 0
    ) -> None:
        self.module = module
        self.match = match
        self.nmatch = nmatch

        if isinstance(self.module, HMMStartEnd):
            self.state = 'm'
        else:
            self.state = None

    def __repr__(self) -> str:
        return str(self.state)

    def determine_state(self, previous: State):
        if self.nmatch == previous.nmatch:
            self.state = 'i'
        else:
            if self.module != PKModule.GAP:
                self.state = 'm'
            else:
                self.state = 'd'


def get_match_states(
    seqs: List[ModuleSequence],
    threshold: Optional[float] = 0.0
) -> List[bool]:
    threshold_count = len(seqs) * threshold
    seqs = [seq._seq for seq in seqs]
    match_states = [
        (sum([module != PKModule.GAP for module, _ in pos]) > threshold_count)
        for pos in list(zip(*seqs))
    ]
    return match_states


def make_logo(path: str, save_to: Optional[str] = None) -> None:
    fasta = parse_fasta(path)
    seqs = [ModuleSequence(record.name, record.seq) for record in fasta]
    match_states = get_match_states(seqs)
    hmm = HMM(match_states)
    hmm.train(seqs)
    hmm.logo(save_to)


def generate_polyketide_backbones(path: str, num: int) -> None:
    fasta = parse_fasta(path)
    seqs = [ModuleSequence(record.name, record.seq) for record in fasta]
    match_states = get_match_states(seqs)
    hmm = HMM(match_states)
    hmm.train(seqs)

    msg = []
    for n in range(1, num + 1):
        module_list = hmm.sample_sequence()
        seq = ModuleSequence(str(n), module_list)
        msg.append(f'>generated_polyketide_backbone_{n}\n{seq}')
    print("\n".join(msg))
