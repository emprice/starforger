const path = require('path');
const webpack = require('webpack');
const TerserPlugin = require('terser-webpack-plugin');

module.exports = {
  mode: 'development',
  devtool: 'source-map',
  entry: path.resolve(__dirname, 'js/index.js'),
  output: {
    path: path.resolve(__dirname, 'public/static'),
    filename: 'starforger.bundle.js',
    library: {
      name: 'starforger',
      type: 'global'
    }
  },
  module: {
    rules: [
      {
        test: /mandelagol.js$/,
        loader: 'exports-loader',
        options: {
          exports: 'default mandelagol'
        }
      },
      {
        test: /mandelagol.wasm$/,
        type: 'asset/resource',
        generator: {
          filename: '[name].wasm'
        }
      }
    ]
  },
  infrastructureLogging: {
    level: 'verbose',
    debug: true
  },
  optimization: {
    minimize: true,
    minimizer: [
      new TerserPlugin({
        minify: TerserPlugin.swcMinify,
        terserOptions: {}
      })
    ]
  },
  resolve: {
    modules: ['node_modules'],
    fallback: {
      'fs': false,
      'path': false
    },
    alias: {
      'autogen': './'
    }
  }
};

// vim: set ft=javascript:
