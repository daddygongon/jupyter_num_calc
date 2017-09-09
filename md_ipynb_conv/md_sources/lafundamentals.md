
線形代数--写像(LAFundamentals)

# 行列と連立方程式

大学の理系で必修なのは微積分と線形代数です．線形代数というと逆行列と固有値の計算がすぐに思い浮かぶでしょう．計算がややこしくてそれだけでいやになります．でも，行列の計算法は一連の手順で記述できるので，Mapleでは微積分とおなじように一個のコマンドで片が付きます．それが3x3以上でも同じです．問題はその意味です．ここでは，線形代数の計算がMapleを使えばどれほど簡単にできるかを示すと共に，線形代数の基本となる概念についてスクリプトと描画を使って，直観的に理解することを目的とします．

先ずは連立方程式から入っていきます．中学の時に

$$
4x = 2
$$
というのを解きますよね．一般的には

$$
\begin {array}{rl}
ax &= b \\
x &= b/a
\end {array}
$$
と書けるというのは皆さんご存知のはず．これと同じようにして連立方程式を書こうというのが逆行列の基本．つまり

$$
\begin {array}{rrl}
2x\, + &5y &=7 \\
4x\, + &y &=5
\end {array}
$$
という連立方程式は，係数から作られる2x2行列を係数行列$A$，左辺の値で作るベクトルを$b$として，

$$
\begin {array}{rll}
Ax &= b & \\
x &= b/A &= A^{-1}b
\end {array}
$$
としたいわけです．

実際にMapleでやってみましょう．行列は英語でMatrixです．
```maple
> restart: A:=Matrix([[2,5],[4,1]]); 
```
$$
A\, := \, \left[ \begin {array}{cc} 2&5\\ 4&1\end {array} \right]
$$
こうして行列を作ります．

```maple
> b:=Vector([7,5]); #(2)ベクトルは英語でVectorです．これで縦ベクトルができます．
```
$$
b\, := \, \left[ \begin {array}{c} 7\\ 5\end {array} \right]
$$
線形代数はlinear algebraと言います．withでLinearAlgebraというライブラリーパッケージを読み込んでおきます．
```maple
> with(LinearAlgebra): 
```
逆行列はmatrix inverseと言います．
```maple
> x0:=MatrixInverse(A).b;
```
行列AのMatrixInverseを求めて，ベクトルbに掛けています．
$$
{\it x0}\, := \, \left[ \begin {array}{c} 1\\ 1\end {array} \right] 
$$
と簡単に求めることができます．


# 掃き出し

係数行列$A$とベクトル$b$を足して作られる行列は拡大係数行列と呼ばれます．Mapleでは，これは
```maple
> <A|b>;
```
$$
\left[ \begin {array}{ccc} 2&5&7\\ 4&1&5\end {array} \right]
$$
として作られます．ここから行列の掃き出し操作をおこなうには，LUDecompositionというコマンドを使います．
```maple
> P,L,U:=LUDecomposition(<A|b>);
```
$$
P,\,L,\,U\, := \, \left[ \begin {array}{cc} 1&0\\ 0&1\end {array} \right] ,\, \left[ \begin {array}{cc} 1&0\\ 2&1\end {array} \right] ,\, \left[ \begin {array}{ccc} 2&5&7\\ 0&-9&-9\end {array} \right]
$$
これは，下三角行列(Lower Triangle Matrix)と上三角行列(Upper Triangle Matrix)に分解(decompose)するコマンドです．$P$行列は置換(permutation)行列を意味します．LUDecompositionだけでは，前進消去が終わっただけの状態です．そこで，後退代入までおこなうには，optionにoutput='R'をつけます．そうすると出力は，
```maple
> LUDecomposition(<A|b>,output='R');
```
$$
\left[ \begin {array}{ccc} 1&0&1\\ 0&1&1\end {array} \right]
$$
で，$b$ベクトルの部分が解になっています．


# 写像

次に，これを2次元上のグラフで見てみましょう．先ず描画に必要なライブラリーパッケージ(plotsおよびplottools)をwithで読み込んでおきます．
```maple
> with(plots):with(plottools): 
```
ベクトルは，位置座標を意味するようにlistへ変換(convert)しておきます．
```maple
> p0:=convert(x0,list); p1:=convert(b,list); 
```
位置p0に円(disk)を半径0.2,赤色で描きます．同じように位置p1に半径0.2，青色でdiskを描きます．もう一つ，p0からp1に向かう矢印(arrow)
を適当な大きさで描きます．後ろの数字をいじると線の幅や矢印の大きさが変わります．
```maple
> point1:=[disk(p0,0.2,color=red), disk(p1,0.2,color=blue)]:
> line1:=arrow(p0,p1,.05,.3,.1 ):
```

これらをまとめて表示(display)します．このとき，表示範囲を-8..8,-8..8とします．
```maple
> display(point1,line1,view=[-2..8,-2..8],gridlines=true);
```

|{{attach_view(LAFundamentalsplot2d1.png,)}}|
|:----|


逆行列は
```maple
> MatrixInverse(A);
```
$$
\left[ \begin {array}{cc} -1/18&5/18\\ 2/9&-1/9\end {array} \right]
%\left[  \begin {array}{cc}  -\frac{1}{18}&   \frac{5}{18}\\
%   \frac{2}{9}&   -\frac{1}{9}\end {array} \right]
$$
で求まります．先ほどの矢印を逆に青から赤へたどる変換になっています．これが，連立方程式を解く様子をグラフで示しています．つまり，行列Aで示される変換によって求まる青点で示したベクトルb(7,5)を指す元の赤点を捜すというものです．答えは(1,1)となります．

では，元の赤点をもう少しいろいろ取って，行列Aでどのような点へ写されるかを見てみましょう．
```maple
> N:=30:point2:=[]:line2:=[]: 
  for k from 0 to N-1 do
    x0:=Vector([sin(2*Pi*k/N),cos(2*Pi*k/N)]); 
    x1:=A.x0; 
    p0:=convert(x0,list);
    p1:=convert(x1,list); 
    point2:=[op(point2),disk(p0,0.05,color=red)];
    point2:=[op(point2),disk(p1,0.05,color=blue)]; 
    line2:=[op(line2),line(p0,p1)];
  end do:
```
N:=30で分割した円周上の点をx0で求めて，point2にその円とそれのＡ.x0を，line2にはその2点を結ぶline(線)を足しています．
使っているコマンドは，先ほどの描画とほぼ同じです．ただし，Mapleスクリプトに特有のidiom(熟語)を使っています．この基本形をとり出すと，
```maple
> list1:=[]; 
  for k from 0 to 2 do 
    list1:=[op(list1),k]; 
  end do; 
  list1;
```
$$
[] \notag \\
[0] \notag \\
[0, 1] \notag \\
[0, 1, 2] \notag \\
[0, 1, 2] \notag
$$
となります．for-loopでkを0から4まで回し，list1に次々と値を追加していくというテクです．

できあがりの次の図を見てください．
```maple
> d:=6: display(point2,line2,view=[-d..d,-d..d]);
```

|{{attach_view(LAFundamentalsplot2d2.png,)}}|
|:----|


何やっているか分かります? 中心の赤点で示される円が，青点で示される楕円へ写されていることが分かるでしょうか．

線形代数の講義で，写像を示すときによく使われるポンチ絵を現実の空間で示すとこのようになります．ポンチ絵では，赤で示した$V$空間が青で示した$W$空間へ行列$A$によって写像され，それぞれの要素$v$が$w$へ移されると意図しています．

|{{attach_view(Projection3-4.png,)}}|
|:----|



# 固有ベクトルの幾何学的意味

では，ここでクイズです．固有ベクトルは上のグラフの何処に対応するか?　ヒントは，
行列Aの固有値，固有ベクトルを$\lambda, x_0$とすると，

$$
A \,x_0 = \lambda \, x_0
$$
が成立する
です．固有値と固有ベクトルはMapleでは以下のコマンドで求まります．
```maple
> lambda,P:=Eigenvectors(A);
```
$$
\lambda,\,P\, := \, \left[ \begin {array}{c} -3\\6\end {array} \right] ,\, \left[ \begin {array}{cc} -1&5/4\\1&1\end {array} \right]
$$
ここではMapleコマンドのEigenvectorsで戻り値を$\lambda$(lambdaと書きます)，$P$に代入しています．この後ろ側にある行列$P$の1列目で構成されるベクトルが固有値-3に対応する固有ベクトル，2列目のベクトルが固有値6に対応する固有ベクトルです．


## 解答 
固有値$\lambda$，固有ベクトル$x_0$の関係式

$$
A \,x_0 = \lambda \, x_0
$$
を言葉で言い直すと，
>固有ベクトル$x_0$は変換行列$A$によって，自分の固有値倍のベクトル$\lambda x_0$に写されるベクトル
となります．つまり変換の図で言うと，
>変換しても方向が変わらない赤点（の方向）
となります．これは図で書くと
```maple
> vv1:=Column(P,1): vv2:=Column(P,2): 
  a1:=vv1[2]/vv1[1]: a2:=vv2[2]/vv2[1]:
  pp1:=plot({a2*x,a1*x},x=-d..d):
```
Columnによって行列の第i列目をとりだし，その比によって直線の傾きを求めています．そうして引いた2本の直線をpp1としてため込んで，先ほど描いた変換の図に加えて表示(display)させます．
```maple
> display(point2,line2,pp1,view=[-d..d,-d..d]);
```

|{{attach_view(LAFundamentalsplot2d3.png,)}}|
|:----|

pp1を入れて描いた直線が引かれた方向ではたしかに変換によっても方向が変わらなさそうに見えるでしょう．

おまけですが，行列の対角化は次のようにしてできます．
```maple
> MatrixInverse(P).A.P;
```
$$
\left[ \begin {array}{cc} -3&0\\ 0&6\end {array} \right]
$$


# 行列式の幾何学的意味

行列Aの行列式($\left|A\right|$あるいはdet$A$と表記)はDeterminantで求まります．
```maple
> Determinant(A);
```
$$
-18
$$
では次のクイズ．先ほど求めた，行列Aの行列式は，どこに対応するでしょう?
以下の(1,0),(0,1)の点を変換した点に原点からベクトルを結んでその意味を説明してください．さらに，そのマイナスの意味は？．
```maple
> point3:=[]:line3:=[]: XX:=Matrix([[1,0],[0,1]]): 
  for i from 1 to 2 do
    x0:=Column(XX,i); x1:=A.x0; 
    p0:=convert(x0,list): 
    p1:=convert(x1,list):
    point3:=[op(point3),disk(p0,0.2,color=red),disk(p1,0.2,color=blue)]:
    line3:=[op(line3),arrow([0,0],p0,.05,.3,.1 ),arrow([0,0],p1,.05,.3,.1 )]:
  end do:
  display(point3,line3,view=[-1..8,-1..8],gridlines=true);
```

|{{attach_view(LAFundamentalsplot2d4.png,)}}|
|:----|



# 行列式が0の写像

では，行列式が０になるというのはどういう状態でしょう? 次のような行列を考えてみましょう．
```maple
> A:=Matrix([[2,1],[4,2]]);
```
$$
$$
この行列式は
```maple
> Determinant(A);
```
$$
0
$$
です．この変換行列で，上と同じように写像の様子を表示させてみましょう．
```maple
> N:=30:point2:=[]:line2:=[]: 
  for k from 0 to N-1 do
    x0:=Vector([sin(2*Pi*k/N),cos(2*Pi*k/N)]); x1:=A.x0; p0:=convert(x0,list);
    p1:=convert(x1,list);
    point2:=[op(point2),disk(p0,0.05,color=red),disk(p1,0.05,color=blue)];
    line2:=[op(line2),line(p0,p1)]; 
  end:
> d:=6: display(point2,line2,view=[-d..d,-d..d]);
```

|{{attach_view(LAFundamentalsplot2d5.png,)}}|
|:----|


わかります？

今回の移動先の青点は直線となっています．つまり，determinantが0ということは，変換すると面積がつぶれるという事を意味しています．平面がひとつ次元を落として線になるということです．

次に，この行列の表わす写像によって原点(0,0)に写される元の座標を求めてみます．連立方程式に戻してみると
```maple
> A.Vector([x,y])=Vector([0,0]);
```
$$
\left[ \begin {array}{c} 2\,x+y\\ 4\,x+2\,y\end {array} \right] = \left[ \begin {array}{c} 0\\ 0\end {array} \right]
$$
となります．とよく見ると，1行目も2行目もおなじ式になっています．2次元正方行列で，行列式が0の時には必ずこういう形になり，直線の式となります．これを表示すると
```maple
> plot([-2*x,-2*x+1,-2*x-1],x=-4..4,y=-4..4);
> plot([2*x],x=-4..4,y=-4..4,color=blue);
```

|{{attach_view(LAFundamentalsKernelImage.png,)}}|
|:----|

左図の赤線となります．この直線上の全ての点が[0,0]へ写されることを確認してください．また，緑の線上の点は全て[1,2]へ写されることが確認できます．
```maple
> A.Vector([-1,2]);
```
$$
\left[ \begin {array}{c} 0\\ 0\end {array} \right]
$$

こうしてすべて調べていけば，左の平面上のすべて点は右の青の直線上へ写されることが分かります．今まで見てきた円と楕円とはまったく違った写像が，行列式が0の行列では起こっていることが分かると思います．右の青線を行列Aによる像(Image, Im$A$と表記)，左の赤線，つまり写像によって[0,0]へ写される集合を核(Kernel, Ker$A$と表記)と呼びます．

これをポンチ絵で描くと，次の通りです．


|像(Image) | 核(Kernel) |
|:----|:----|
|{{attach_view(Projection1.png,LAFundamentals)}}|{{attach_view(Projection2.png,LAFundamentals)}}|




# 全単射

行列$A$による写像を$f$として，赤点に限らず元の点の集合を$V$, 移った先の点の集合を$W$とすると，

$$
f: V \rightarrow W
$$
と表記されます．$v,w$を$V,W$の要素としたとき，異なる$v$が異なる$w$に写されることを単射，全ての$w$に対応する$v$がある写像を全射と言います．全単射，つまり全射でかつ単射，だと要素は一対一に対応します．先ほどのAは全射でもなく，単射でもない例です．

行列式が0の場合の写像は単射ではありません．このとき，逆写像が作れそうにありません．これを連立方程式に戻して考えましょう．もともと，

$$
v = A^{-1} w
$$
の解$v$は点$w$が写像$A$によってどこから写されてきたかという意味を持ちます．逆写像が作れない場合は，連立方程式の解はパラメータをひとつ持った複数の解(直線)となります．これが係数行列の行列式が0の場合に，連立方程式の解が不定となる，あるいは像がつぶれるという関係です．

行列の次元が高い場合には，いろいろなつぶれかたをします．行列の階数と次元は
```maple
> Rank(A); 
  Dimension(A);
```
$$
1 \notag \\
2, 2 \notag
$$
で求まります．

Aをm行n列の行列とするとき，
>Rank(*A'') = Dimension (Im ''A*)

>Dimension (Ker *A'') = ''n'' - Rank(''A*) 

が成立し，これを次元定理といいます．
全射と単射の関係は，下の表のような一変数の方程式での解の性質の拡張と捉えることができます．

### caption:代数方程式$a x =b$の解の存在性．

|一意|0$ ||
|:----|:----|:----|:----|
|不定|$a=0, b=0$ |解は無数 |
|不能|0$ ||




|m x n行列A|全射でない(Im $A < m$), 値域上にあるときのみ解が存在|全射(Im $A =m$), 解は必ず存在|
|:----|:----|:----|:----|
| 0$), 解は複数 ||{{attach_view(Projection3-2.png,LAFundamentals)}}|
|単射(Ker $A = 0$), 解はひとつ|{{attach_view(Projection3-3.png,LAFundamentals)}}|{{attach_view(Projection3-4.png,LAFundamentals)}}|





# 課題

1.  下の図は
$A\, := \, \left[ \begin {array}{cc} 3&2/3\\ 2/3&2\end {array} \right]$を用いて変換される像を表わしている．この行列の固有値，行列式が何処に対応するか説明せよ．また，固有ベクトルの方向を記せ．(2007年度期末試験)

|{{attach_view(LAFundamentalsplot2d8.png,)}}|
|:----|
