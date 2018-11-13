数式処理テキスト
================

基礎
----

-   [first\_leaf.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/first_leaf.ipynb)
-   [functions.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/functions.ipynb)
-   [equals.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/equals.ipynb)
-   課題
    -   [group\_works\_1.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_1.ipynb)
    -   [group\_works\_1\_ans.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_1_ans.ipynb)

微積
----

-   [differential.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/differential.ipynb)
-   [integral.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/integral.ipynb)
-   課題
    -   [group\_works\_2.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_2.ipynb)
    -   [group\_works\_2\_ans.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_2_ans.ipynb)

線形代数
--------

-   [linear\_algebra\_scipy.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/linear_algebra_scipy.ipynb)
-   [linear\_algebra\_sympy.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/linear_algebra_sympy.ipynb)
-   課題
    -   [file:group\_works\_3.pdf](group_works_3.pdf)
    -   [group\_works\_3\_breast\_cancer.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_3_breast_cancer.ipynb)
    -   [group\_works\_3\_ans.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_3_ans.ipynb)
        -   [file:validate\_A.data](validate_A.data),[file:validate\_b.data](validate_b.data),[file:test\_A.data](test_A.data),[file:test\_b.data](test_b.data)
    -   [file:group\_works\_4.pdf](group_works_4.pdf)

式変形
------

-   [equation\_manipulation.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/equation_manipulation.ipynb)
-   課題
    -   [group\_works\_5.pdf](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_5.pdf)
    -   [group\_works\_5\_ans.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/group_works_5_ans.ipynb)

その他
------

-   [python\_files\_and\_dir.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/python_files_and_dir.ipynb)
-   [miscellaneous.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/miscellaneous.ipynb)
-   [cg.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/cg.ipynb)
-   [linear\_solve.ipynb](https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/linear_solve.ipynb)

jupyter notebook 授業でのQ & A
==============================

授業プリントをnbviewerで開けるとnot foundがでる
:   browserを再起動

家でやってきたのをどうすれば見れる
:   z:/に入れる

前に同じマシンでやったのは?
:   c:/users/自分のID

人シスとか，4階でやったのは？
:   残念ながら，あっちのはlogoutしたら消えます．

課題のコメントは?
:   cellの属性をmarkdownに変えて入力，latex記法が使えます．

jupyter notebookが反応しない
:   kernel, interruptかrestartを掛ける

nbviewerからipynb
:   down loadできます

nbviewerの表示が古い
:   原因はわかりません．しばらくしてから再度アクセスしてみてください．

jupyter notebookのinit時の問題点
================================

起動時のhome directory
----------------------

今はz:/です．先週まではlocal:/\[USER ID\]でした．
IV,V号館の演習室のPCでは4月いっぱいはlocal diskになっているので，
logoutすると消えちゃいます．USBに写しておいてください．来月から同じになります．

起動
----

書類のダブルクリックで起動できない．

printout
--------

IEだけでなく，Chromeでも後ろが切れるというバグがありました．
起動するbrowserをFireFoxにしたら大丈夫です．

1.  Jupyter notebookのmenu barのFile-&gt;Print Previewで作って，
2.  FireFoxの右上のツールボタンから印刷を選んでください．

集約して出すのを忘れずに．

集約は，出てきた画面の左上の印刷ボタンを選んで，「プロパティ」から
ワンクリック設定アイコンで集約を選択して印刷．

-   出力が半分切れるような表示になる．

init\_session()でMultipleInstanceErrorがでる :: Googleだと
----------------------------------------------------------

<https://github.com/sympy/sympy/issues/13319> とのことだが...
三月ごろ修正パッチが貼られたversionがrelease?,
init\_print()で代用，あとkernel restart

init\_session()
---------------

-   init\_session() =&gt; init\_printing()

すると綺麗に出力してくれる。

-   出力が混乱してきたときは、Kernel restartで初期化するべし。

解けない問題リスト
------------------

-   16\_pre 1(a) : (x\*\*3+1)\*\*(1/3) sqrt-&gt;cbrt (5523),
-   16\_pre 2(a) : int(1/(cos(x)\*\*2 + 4\*sin(x)\*\*2),x)
-   13-2 1(a) : int(1/(2+cos(t)),t)
-   16-1-pair 2(b)
-   rational
-   aspect ratio = 1
-   log plot, sympyはOKだが，numpyだと悲惨．

