
微積分

# 微分(Diff)-I


## 単純な微分(diff) 
単純な一変数関数の一次微分は，以下の通り．
```maple
> diff(x^2-3*x+2,x);	#res: 2x-3
```
高次の微分は，微分変数を必要なだけ並べる．
```maple
> diff(sin(x),x,x);	#res: -sin(x)
```
さらに高次では次のように\$を使った記法が便利．これはxについての3次微分を表わす．
```maple
> diff(x^4,x$3);	#res: 24x
```


#### 偏微分(PartialDiff) 
複数の変数を持つ多変数の関数では，微分する変数を明示すれば偏微分が求められる．
```maple
> eq1:=(x+y)/(x*y);
> diff(eq1,x);
```
$$
eq1 := \,{\frac {x+y}{xy}}
$$
$$
{\frac {1}{xy}}-{\frac {x+y}{{x}^{2}y}}
$$



## 例題:関数の微分と増減表 

次の関数とその1次導関数を同時にプロットし概形を確認し，さらに増減表を求めよ．


$$
\frac {x}{{x}^{2}-2x+4}
$$


#### 解答例 

```maple
> f0:=unapply(x/(x^2-2*x+4),x):
> df:=unapply(diff(f0(x),x),x);
> plot([f0(x),df(x)],x);
```

$$
{\it df}\, := \,x\mapsto  \left( {x}^{2}-2\,x+4 \right) ^{-1}-{\frac {x \left( 2\,x-2 \right) }{ \left( {x}^{2}-2\,x+4 \right) ^{2}}}
$$


|{{attach_view(Diffplot2d1.png,)}}|
|:----|




## 例題:接線(Tangent) 
次の関数の$ x=3$での接線を求め，2つの関数を同時にプロットせよ．

$$
y={x}^{3}-2\,{x}^{2}-35\,x
$$


#### 解答例 
与関数をf0と定義．
```maple
> f0:=unapply(x^3 - 2*x^2 - 35*x,x);
```
$$
{\it f0}\, := \,x\mapsto {x}^{3}-2\,{x}^{2}-35\,x
$$
微分関数をdfと定義
```maple
> df:=unapply(diff(f0(x),x),x);
```
$$
{\it df}\, := \,x\mapsto 3\,{x}^{2}-4\,x-35
$$
接点(x0,f0(x0))で傾きdf(x0)の直線をf1と定義．
```maple
> x0:=3;
> eq1:=df(x0)*(x-x0)+f0(x0);
> f1:=unapply(eq1,x);
```
$$
{\it x0}\, := \,3 \notag \\
{\it eq1}\, := \,-20\,x-36 \notag \\
{\it f1}\, := \,x\mapsto -20\,x-36 \notag
$$
2つの関数を同時にプロット．
```maple
> plot([f0(x),f1(x)],x=-5..5);
```

|{{attach_view(Diffplot2d2.png,)}}|
|:----|



# 微分(Diff)-II


## 級数展開(series) 
Taylor級数は以下のようにして，中心点(x=a)，次数(4次)を指定する．
```maple
> t1:=series(sin(x),x=a,4);
```
$$
t1 := \sin(a)+\cos(a)(x-a)-\frac{1}{2}\sin(a)(x-a)^2-\frac{1}{6}\cos(a)(x-a)^3+O((x-a)^4)
$$
```maple
> e1:=convert(t1,polynom);
> f1:=unapply(e1,x);
```
$$
t1 := \sin(a)+\cos(a)(x-a)-\frac{1}{2}\sin(a)(x-a)^2-\frac{1}{6}\cos(a)(x-a)^3 \notag \\
{\it f1}\, := \,x\mapsto 
\sin(a)+\cos(a)(x-a)-\frac{1}{2}\sin(a)(x-a)^2-\frac{1}{6}\cos(a)(x-a)^3 \notag
$$



## 全微分(D) 
全微分を計算するときは，Dを用いる．
```maple
> f:=unapply(x^4*exp(-y^2),(x,y));
> D(f(x,y));
> (D@@2)(f(x,y));
```
$$
f\, := \,( {x,y} )\mapsto {x}^{4}\exp(-{y}^{2}) \notag \\
4\, {D} \left( x \right) {x}^{3}\exp(-{y}^{2})+{x}^{4} {D} \left( \exp(-{y}^{2}) \right) \notag \\
4\, \left( D^{ \left( 2 \right) } \right)  \left( x \right) {x}^{3}\exp(-{y}^{2})+12\, \left(  {D} \left( x \right)  \right) ^{2}{x}^{2}\exp(-{y}^{2})+8\, {D} \left( x \right) {x}^{3} {D} \left( \exp(-{y}^{2}) \right) +{x}^{4} \left( D^{ \left( 2 \right) } \right)  \left( \exp(-{y}^{2}) \right) \notag
$$

ここで，D(x)などはxの全微分を表わす．これは，x,yを変数としているので
```maple
> diff(x,x);
> diff(exp(-y^2),y);
```
$$
1 \notag \\
-2\,y\exp(-{y}^{2}) \notag
$$
であるがMapleには分からない．そこで全微分の最終形を得るには，あらかじめD(x)などの結果を求めておき，subsで明示的に代入する必要がある．
```maple
> dd:=D(f(x,y)):
> eqs:={D(x)=diff(x,x),D(exp(-y^2))=diff(exp(-y^2),y)};
> subs(eqs,dd);
```
$$
{\it eqs}\, := \, \left\{  {D} \left( x \right) =1, {D} \left( \exp(-{y}^{2}) \right) =-2\,y\exp(-{y}^{2}) \right\}  \notag \\
4\,{x}^{3}\exp(-{y}^{2})-2\,{x}^{4}y\exp(-{y}^{2}) \notag
$$



## 複合関数の微分 
```maple
> diff(f(x)*g(x),x);
> diff(f(g(x)),x);
```
$$
\left( {\frac {d}{dx}}f \left( x \right)  \right) g \left( x \right) +f \left( x \right) {\frac {d}{dx}}g \left( x \right) \notag \\
\mbox {D} \left( f \right)  \left( g \left( x \right)  \right) {\frac {d}{dx}}g \left( x \right) \notag
$$

```maple
> f:=x->exp(x);
> g:=x->cos(x);
> diff(f(x)*g(x),x);
> diff(f(g(x)),x);
```
$$
f\, := \,x\mapsto \exp(x) \notag \\
g\, := \,x\mapsto \cos \left( x \right)  \notag \\
\exp(x)\cos \left( x \right) -\exp(x)\sin \left( x \right)  \notag \\
-\sin \left( x \right) \exp(\cos x ) \notag
$$



### 微分に関する課題 
1. 
次の関数を微分せよ．

i)$ {x} \log x$, 
ii) $ \frac{1}{  \left( 1+x \right) ^{3}}$, 
iii) $ \sqrt{4\,x+3}$, 
iv) $ \frac{1}{ a^2+ \left( x-x_0 \right)^2 }$

1. 
次の関数の1次から5次導関数を求めよ．

i) $\sin^2 x$, 
ii) ${e}^{x}$

1.  以下の関数をx0まわりで３次までテイラー展開し，得られた関数ともとの関数をプロットせよ．さらに5次まで展開した場合はどう変化するか．

i) $ y=\sin x, x_0=0 $, 
ii) $ y=\cos x, x_0=\frac{\pi}{2}$

1. 
(発展課題）$f \left( x,y \right) ={e}^{x}{\it log} \left( 1+y \right) $
を$ x=0,\,y=0$のまわりで3次まで展開せよ．




### Diff 
1.  
```maple
> diff(x*log(x),x);
> diff(1/(1 + x)^3,x);
> diff(sqrt(4*x + 3),x);
>diff(1/(a^2+(x-x0)^2),x);
```
$$
\ln  \left( x \right) +1 \notag \\
-3\, \left( 1+x \right) ^{-4} \notag \\
2\, \left(  \sqrt{4\,x+3} \right) ^{-1} \notag \\
-{\frac {2\,x-2\,{\it x0}}{ \left( {a}^{2}+ \left( x-{\it x0} \right) ^{2} \right) ^{2}}} \notag
$$

1. 
```maple
> diff(sin(x)^2,x);
> diff(sin(x)^2,x$2);
```

$$
2\,\sin \left( x \right) \cos \left( x \right) \notag \\
2\, \left( \cos \left( x \right)  \right) ^{2}-2\, \left( \sin \left( x \right)  \right) ^{2}\notag
$$
```maple
以下略
```


1. 
先ず与関数をf0と定義
```maple
> f0:=unapply(sin(x),x);
```
$$
{\it f0}\, := \,x\mapsto \sin \left( x \right)
$$
テイラー展開した結果をeq1とする．関数として定義するためにeq1を多項式に変換し(convert)，unapplyをかける．
```maple
> eq1:=series(f0(x),x=0,3);
> f1:=unapply(convert(eq1,polynom),x);
```
$$
{\it eq1}\, := \,x+O \left( {x}^{3} \right) \notag \\
{\it f1}\, := \,x\mapsto x \notag
$$
５次についても同様
```maple
> eq2:=series(f0(x),x=0,5);
> f2:=unapply(convert(eq2,polynom),x);
```
$$
{\it eq2}\, := \,x-1/6\,{x}^{3}+O \left( {x}^{5} \right) \notag \\
{\it f2}\, := \,x\mapsto x-1/6\,{x}^{3} \notag
$$
3つの関数を同時プロット
```maple
> plot([f0(x),f1(x),f2(x)],x=-Pi..Pi);
```

|{{attach_view(Diffplot2d3.png,)}}|
|:----|


```maple
> series(f0(x),x=Pi/2,3)
```
以外は前問とおなじ．

1. 
```maple
f:=unapply(exp(x)*log(1+y),(x,y));
```
$$
f\, := \,( {x,y} )\mapsto \exp(x)\ln  \left( 1+y \right)
$$
```maple
eq1:=series(series(f(x,y),x=0,3),y=0,3);
```
$$
{\it eq1}\, := \,O \left( {x}^{3} \right) + \left( 1+1/2\,{x}^{2}+x \right) y+ \left( -1/2\,x-1/2-1/4\,{x}^{2} \right) {y}^{2}\\
\mbox{}+O \left( {y}^{3} \right)
$$

```maple
> g:=unapply(convert(convert(eq1,polynom),polynom),(x,y));
```
$$
g\, := \,( {x,y} )\mapsto  \left( 1+1/2\,{x}^{2}+x \right) y+ \left( -1/2\,x-1/2-1/4\,{x}^{2} \right) {y}^{2}
$$
```maple
> plot3d([f(x,y),g(x,y)],x=-1..1,y=-1..1,axes=box);
```

|{{attach_view(Diffplot3d4.png,)}}|
|:----|




### きれいな表示 

以下のようにすると表示がきれい．
```maple
> f:=unapply(x^4*exp(-y^2),(x,y));
> d:=Diff(f(x,y),x);
> d=value(d);
```
$$
f\, := \,( {x,y} )\mapsto {x}^{4}\exp(-{y}^{2}) \notag \\
d\, := \,{\frac {\partial }{\partial x}} \left( {x}^{4}\exp(-{y}^{2}) \right) \notag \\
{\frac {\partial }{\partial x}} \left( {x}^{4}\exp(-{y}^{2}) \right) =4\,{x}^{3}\exp(-{y}^{2})\notag
$$
