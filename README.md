The multi-reference density functional theory with self consistent field (MR-DFT-SCF) method
                              (Version 1.0 alpha)


 Main auther       :  Zexing  Qu  (Jilin Univ.)            zxqu@jlu.edu.cn

 Main contributer  :  Yingjin Ma  (CNIC@CAS)               yingjin.ma@sccas.cn
                      Adam        (Minnesota Univ.)


Command for git:

git status
git checkout -- [file name]
git branch
git checkout [remote branch]
git pull

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
