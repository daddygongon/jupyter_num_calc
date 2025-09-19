## VS Code + Jupyter Notebook (ipynb) セットアップ手順（macOS, Homebrew Python, venv利用）

## 1. 仮想環境（venv）の作成

ターミナルで作業ディレクトリに移動し、仮想環境を作成します。

```sh
cd /Users/bob/Desktop/Lectures/math_python/jupyter_num_calc/symbolic_math
python3.12 -m venv venv
```

## 2. 仮想環境をアクティブ化

```sh
source venv/bin/activate.fish
```

## 3. 必要なパッケージのインストール

```sh
pip install --upgrade pip
pip install ipykernel sympy matplotlib numpy
```

## 4. VS Codeでノートブックを開く

- VS Codeで `.ipynb` ファイルを開きます。

## 5. カーネル（インタープリタ）の選択

- ノートブック右上の「カーネル」または「Python 3.12.11」などと表示されている部分をクリック。
- `venv`（例: `/Users/bob/Desktop/Lectures/math_python/jupyter_num_calc/symbolic_math/venv/bin/python`）を選択。

  ※ 表示されない場合は、コマンドパレット（`Cmd+Shift+P`）で  
  `Notebook: Select Notebook Kernel` を実行し、`venv` を選択。

## 6. ノートブックを実行

- セルを実行して、エラーが出ないことを確認します。

---

## トラブルシューティング
- HomebrewのPythonでは `pip install --user` などが制限されているため、**必ず仮想環境を使う**こと。
- `ipykernel` のエラーが出る場合は、カーネル（インタープリタ）が仮想環境のものになっているか確認。
- `import sys; print(sys.executable)` で現在のPythonパスを確認できる。

---

## 参考

- VS Code拡張機能「Python」「Jupyter」がインストールされていることを確認してください。