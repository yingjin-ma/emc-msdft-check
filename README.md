The multi-reference density functional theory with self consistent field (MR-DFT-SCF) method
                              (Version 1.0 alpha)
                                                Zexing Qu (Jilin University)
                                                Email: zxqu@jlu.edu.cn
                                                                   2017.9.27
Command for git:

CASE 1:
Git global setup

git config --global user.name "Zexing Qu"
git config --global user.email "zxqu@jlu.edu.cn"

CASE 2:
Create a new repository

git clone http://zxqu@git.vlcc.cn/zxqu/MRDFTSCF.git
cd MRDFTSCF
touch README.md
git add README.md
git commit -m "add README"
git push -u origin master

CASE 3:
Existing folder

cd existing_folder
git init
git remote add origin http://zxqu@git.vlcc.cn/zxqu/MRDFTSCF.git
git add .
git commit
git push -u origin master

CASE 4:
Existing Git repository

cd existing_repo
git remote add origin http://zxqu@git.vlcc.cn/zxqu/MRDFTSCF.git
git push -u origin --all
git push -u origin --tags
