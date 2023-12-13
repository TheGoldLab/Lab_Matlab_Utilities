from typing import Any

import numpy as np

from pyramid.trials.trials import Trial, TrialCollecter


class TACPCollecter(TrialCollecter):
    """A simple enhancer that computes trials after change-point for each trial"""

    def __init__(self) -> None:
        self.tacp = 0

    def __eq__(self, other: object) -> bool:
        """Compare by attribute, to support use of this class in tests."""
        if isinstance(other, self.__class__):
            return self.max_start_time == other.max_start_time
        else:  # pragma: no cover
            return False

    def __hash__(self) -> int:
        """Hash by attribute, to support use of this class in tests."""
        return self.max_start_time.__hash__()

    def collect(
        self,
        trial: Trial,
        trial_number: int,
        experiment_info: dict,
        subject_info: dict
    ) -> None:
        # Store the max start time over all trials.
        self.max_start_time = max(self.max_start_time, trial.start_time)

    def enhance(
        self,
        trial: Trial,
        trial_number: int,
        experiment_info: dict[str: Any],
        subject_info: dict[str: Any]
    ) -> None:
        # Compute the start time of this trial as a percentage of the whole session.
        percent_complete = 100 * trial.start_time / self.max_start_time
        trial.add_enhancement("percent_complete", percent_complete, "value")
