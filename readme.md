#Instructions

1. ```cd static```
2. ```npm init```
3. ```npm i webpack --save-dec```
4. 
```
const webpack = require('webpack');

const config = {
    entry: __dirname + '/js/index.jsx', 
    output: {
        path: __dirname + '/dist',
        filename: 'bundle.js',

    },
    resolve: {
        extensions: ['.js', '.jsx', '.css']
    },
    module: {
        rules: [
          {
            test: /\.jsx?/,
            exclude: /node_modules/,
            use: 'babel-loader'
          }
        ]
      },
}

module.exports = config;
```
5. 
``` 
"scripts": {
    "build": "webpack -p --progress --config webpack.config.js",
    "dev-build": "webpack --progress -d --config webpack.config.js",
    "test": "echo \"Error: no test specified\" && exit 1",
    "watch": "webpack --progress -d --config webpack.config.js --watch"
  },
``` 
6. 
``` 
npm i babel-core babel-loader babel-preset-es2015 babel-preset-react --save-dev
npm i @babel/preset-env @babel/preset-react --save-dev
```
7. See package json
8. See index.html
9. Create index.jsx file
10. ```npm run watch```
11. ```npm i react react-dom --save-dev```
12. cd to app.py
13. ```python3 -m venv env```
14. ```source env/bin/activate ```
15. ``` pip install flask ```
16. See app.py
17. ```python app.py```
