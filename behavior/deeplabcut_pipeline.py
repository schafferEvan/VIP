#implementation of the DeepLabCut Model
#This code is a direct translation of the Jupyter Notebook example labelled Demo_yourowndata.ipnb from DeepLabCut examples.
#import sys
#sys.path.append('/mnt/home/evanschaffer/anaconda3/envs/deeplabcut')
import os
os.environ['DLClight'] = 'True'
import deeplabcut


def create_directory(task,experimenter,video,working_directory):
    path=deeplabcut.create_new_project(task,experimenter,video, working_directory,copy_videos=True) #change the working directory to where you want the folders created.
    return path

def extract_frames(path_config_file):
    # %matplotlib inline
    deeplabcut.extract_frames(path_config_file,'automatic','uniform',crop=True) #there are other ways to grab frames, such as by clustering 'kmeans'; please see the paper.

def label_frames(path_config_file):
    # %gui wx
    deeplabcut.label_frames(path_config_file)

    deeplabcut.check_labels(path_config_file) #this creates a subdirectory with the frames + your labels
def create_training_set(path_config_file):
    deeplabcut.create_training_dataset(path_config_file)
def train_network(path_config_file):
    deeplabcut.train_network(path_config_file)
def evaluate_network(path_config_file):
    deeplabcut.evaluate_network(path_config_file)
# path_config_file = 'test2/Walking-Benjamin-2019-05-31/config.yaml'

def analyze_video(path_config_file,videofile_path):
    deeplabcut.analyze_videos(path_config_file,videofile_path)
def extract_outlier(path_config_file):
    deeplabcut.extract_outlier_frames(path_config_file,['/videos/video3.avi'])
def refine_labels(path_config_file):
    # %gui wx
    deeplabcut.refine_labels(path_config_file)
    #Once all folders are relabeled, check them and advance. See how to check labels, above!
    deeplabcut.merge_datasets(path_config_file)
    create_training_set(path_config_file)

def create_labeled_video(path_config_file,videofile_path):
    deeplabcut.create_labeled_video(path_config_file,videofile_path)
def plot_trajectories(path_config_file,videofile_path):
    # %matplotlib notebook #for making interactive plots.
    deeplabcut.plot_trajectories(path_config_file,videofile_path)


# MAIN FUNCTION
"""
stage 1 sets up frames and labels.
    - must run with pythonw
    - must have the hui active
stage 2 trains the network and makes labeled video with graphs.
    - can run with regular python
"""
def main():
    stage = 2
    # print('\n'+100*"-"+'\n')
    # import sys
    # sys.path.append('C:/Users/Benjamin/Desktop/DeepLabCut/deeplabcut')
    # import deeplabcut
    # import 
    task='Walking' # Enter the name of your experiment Task
    experimenter='Benjamin' # Enter the name of the experimenter
    video=['/mnt/ceph/users/evanschaffer/data/deeplabcut/cluster_t1/train_videos/fly1_run2.mp4','/mnt/ceph/users/evanschaffer/data/deeplabcut/cluster_t1/train_videos/fly2_run2.mp4'] # Enter the paths of your videos you want to grab frames from.
    working_directory='/mnt/ceph/users/evanschaffer/data/deeplabcut/cluster_t1'
    videofile_path = ['/mnt/ceph/users/evanschaffer/data/deeplabcut/cluster_t1/test_videos/fly3_run2.mp4']

    # cpath = input("Create new path? (y/n): ")
    # if cpath == 'y':
    #     path_config_file = create_directory(task,experimenter,video,working_directory)
    # else:
    path_config_file = '/mnt/ceph/users/evanschaffer/data/deeplabcut/cluster_t1/Walking-Benjamin-2019-06-04/config.yaml'
    # proceed = input('\nSet up your config.yaml file with the correct labels: \n')

    #frame setup
    if stage == 1:
        import matplotlib
        matplotlib.use('Agg')
        extract_frames(path_config_file)
        label_frames(path_config_file)
        proceed = input('\nPress key to build network: \n')


    else:
        # start network building
        # create_training_set(path_config_file) 
        # proceed = input('\nTime to set of your train .yaml file. Press any key when ready: \n') 
        
        # start training
        train_network(path_config_file)
        evaluate_network(path_config_file)
        proceed = input('\nPress key to create labeled video with trajectories: \n')

        # #start video analysis
        # analyze_video(path_config_file,videofile_path)
        # create_labeled_video(path_config_file,videofile_path)
        # plot_trajectories(path_config_file,videofile_path)

if __name__ == '__main__':
    main()


