import os
import sys
import pickle


# Set memory
os.environ['MALLET_MEMORY'] = '200G'

from pycisTopic.lda_models import run_cgs_models_mallet

# Load cistopic object
out_dir = "outs"
with open(os.path.join(out_dir, "cistopic_obj.pkl"), "rb") as f:
    cistopic_obj = pickle.load(f)

print("✓ Loaded cistopic object")

# Configure Mallet
mallet_path = "/projectnb/paxlab/presh/software/Mallet-202108/bin/mallet"

# Create temp directory with write permissions
tmp_path = os.path.join(out_dir, "mallet_tmp")
save_path = os.path.join(out_dir, "lda_models")
os.makedirs(tmp_path, exist_ok=True)
os.makedirs(save_path, exist_ok=True)

print("Starting LDA model training...")

try:
    models = run_cgs_models_mallet(
        cistopic_obj,
        n_topics=[10, 15, 20, 25, 30, 35, 40, 45, 50],
        n_cpu=12,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path=tmp_path,  # Use project temp
        save_path=save_path,
        mallet_path=mallet_path,
    )
    
    print("✓ Models trained successfully")
    
    # Save models
    with open(os.path.join(out_dir, "lda_models2.pkl"), "wb") as f:
        pickle.dump(models, f)
    
    print("✓ Models saved")
    
except Exception as e:
    print(f"Error during model training: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)