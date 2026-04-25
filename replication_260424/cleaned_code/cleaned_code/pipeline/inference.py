"""ST3DCN model loading + batch inference glue.

The paper's repo must be on sys.path (for `ST3DCN_Utils` import) and
CUDA_VISIBLE_DEVICES + LD_LIBRARY_PATH must be set before calling `load_model`.

Typical usage:
    os.environ['CUDA_VISIBLE_DEVICES'] = '2'
    os.environ['LD_LIBRARY_PATH'] = '/usr/local/cuda-11.2/lib64:' + os.environ.get('LD_LIBRARY_PATH','')
    sys.path.insert(0, '/path/to/ST3DCN')
    model = load_model('/path/to/ST3DCN_Model.h5')
    prob = predict(model, crop)
"""
import numpy as np


THRESHOLD_HCC_BINARY = 0.8  # Paper spec (binary HCC vs not-HCC)


def load_model(weights_path: str, width: int = 70, height: int = 70, depth: int = 70,
                batch_size: int = 16, factor: int = 8, num_class: int = 2):
    """Load ST3DCN with pretrained weights."""
    import tensorflow as tf
    for g in tf.config.list_physical_devices('GPU'):
        tf.config.experimental.set_memory_growth(g, True)
    from ST3DCN_Utils import multi_scale_get_model_DCN
    model = multi_scale_get_model_DCN(width=width, height=height, depth=depth,
                                        batch_size=batch_size, factor=factor,
                                        num_class=num_class)
    model.load_weights(weights_path)
    return model


def predict(model, crop: np.ndarray) -> float:
    """Run single-lesion inference. `crop` shape (D, H, W) float32 [0,1].
    Model expects (batch, D, H, W, 1); returns scalar p(HCC)."""
    x = crop[None, :, :, :, None]
    return float(model.predict(x, verbose=0).ravel()[0])


def binary_pred(prob: float, threshold: float = THRESHOLD_HCC_BINARY) -> int:
    return int(prob >= threshold)
