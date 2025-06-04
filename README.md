# F90_jacobianEigenvalue

FK3の計算結果に対してCEMA (Chemically Explosive Mode Analysis)を適用するためのコード．  
Paraview用の可視化スクリプトpvをもとに作成しているので，`$PROJECT_DIR/data/F90_jacobianEigenvalue`のように配置して，vtsファイルを作成可能．  

## nh3_reduced branch
簡略化したTamaokiメカニズムを使用する

## 使い方

解析的なヤコビアンを取得するために[pyJach](https://github.com/SLACKHA/pyJac)で生成されたコードを`src/c`に配置して使用している．そのため,初めにpyJacのサイトの指示に従って反応機構をcのコードに変換する処理が必要．これによって固有値と爆発指数(Explosion Index: EI)が計算可能になる．

参与指数(Participation Index: PI)を取得するためには，上記の処理に加えて`ipynb/get_stoichiometric_coefficients.ipynb`を用いて，反応速度に関するパラメータをの前処理を行う必要がある．

実行は以下のコマンドにより行う．

```bash
make
bin/main
```

実行後にxy_planeにvtsファイルが出力されるのでこれをParaviewで読み込むことができる．

EIとPIのラベルを付与するためには，`ipynb/annotate_vts.ipynb`を実行すると_annotation.vtsというファイルがxy_planeに出力される．