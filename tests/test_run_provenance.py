"""The training script records the commit it runs from, so a result names the code that produced it."""
import importlib.util
import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_training_script():
    specification = importlib.util.spec_from_file_location("run_main_model", REPOSITORY_ROOT / "experiments" / "run_main_model.py")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_git_provenance_names_the_head_commit_and_lists_tracked_changes() -> None:
    provenance = load_training_script().git_provenance()
    assert re.fullmatch(r"[0-9a-f]{40}", provenance["commit"])
    assert isinstance(provenance["tracked_changes"], list)


def test_optimizer_groups_give_links_and_leaks_their_own_time_scales() -> None:
    from argparse import Namespace

    import torch

    from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead

    training_script = load_training_script()
    encoder = torch.nn.Linear(2, 2)
    head = NoisyOrPathwayModuleHead(num_graph_nodes=4, node_state_dim=2, num_pathway_modules=2, num_symptoms=3)
    groups = training_script.optimizer_parameter_groups(encoder, head, Namespace(link_learning_rate=0.02, leak_learning_rate=0.0002, module_bias_learning_rate=0.0, gate_learning_rate=0.05))
    rates_by_parameter = {id(parameter): group.get("lr") for group in groups for parameter in group["params"]}
    assert rates_by_parameter[id(head.module_symptom_link_logit)] == 0.02
    assert rates_by_parameter[id(head.symptom_leak_logit)] == 0.0002
    assert rates_by_parameter[id(head.support_gate.log_alpha)] == 0.05
    assert rates_by_parameter[id(head.module_readout_bias)] is None  # main learning rate
    assert len(rates_by_parameter) == len(list(encoder.parameters())) + len(list(head.parameters()))
    single_group = training_script.optimizer_parameter_groups(encoder, head, Namespace(link_learning_rate=0.0, leak_learning_rate=0.0, module_bias_learning_rate=0.0))
    assert len(single_group) == 1
