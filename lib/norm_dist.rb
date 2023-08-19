# frozen_string_literal: true

require_relative "norm_dist/version"
require "distribution"

module NormDist
  class Error < StandardError; end
  
  # 自力で正規分布関数を構築
  def normal_pdf_calc(x, mu, sigma)
    expo = -(((x - mu) ** 2) / (2 * (sigma ** 2)))
    return 1.0 / Math::sqrt(2.0 * Math::PI * (sigma ** 2.0)) * Math::exp(expo)
  end

  # 自力で累積正規分布関数を構築
  # calc_range_multiple: 積分開始点の平均値からの距離(σのcalc_range_multiple倍として定義)
  # steps: 積分ステップ数
  def normal_cdf_calc(x, mu, sigma, calc_range_multiple: 10.0, steps: 25000)
    ### 積算する範囲設定 この調整がむずい
    len = sigma*calc_range_multiple < (x-mu).abs * 2 ? (x-mu).abs * 2 : sigma * calc_range_multiple
    step_len = (x-(mu-len)) / steps
    calc_range = ((mu-len)..x).step(step_len)
    ###
    
    cdf = calc_range.inject(0.0) do |result, cx|
      pdf1 = normal_pdf(cx-step_len, mu, sigma)
      pdf2 = normal_pdf(cx, mu, sigma)
      plus = (pdf1 + pdf2) * step_len / 2.0
      result += plus
    end
    return cdf
  end

  ## 標準化を使用するパターン
  def normal_pdf(x, mu, sigma)
    Distribution::Normal.pdf((x - mu) / sigma) / sigma
  end

  ## 多分、正規分布を標準化とか軒並み言ってるのは累積分布のことなんだろう。
  ## なにしろ標準分布表の値のことなんだから。
  def normal_cdf(x, mu, sigma)
    Distribution::Normal.cdf((x - mu) / sigma)
  end

  ## 平均値の区間推定
  def estimate(df, colname, n: 20, print_data: true)
    # サンプリング
    samples = df.row[*(df.index.to_a.sample(n))]

    # 定数の計算
    m = samples[colname].mean
    v = samples[colname].variance
    sig = v ** 0.5

    # 信頼区間95%の推定
    area =  [m - 1.96 * sig / Math::sqrt(n), m + 1.96 * sig / Math::sqrt(n)]
    print_data && print(area)
    
    # 検証
    area_range = area[0]..area[1]
    true_ave = df[colname].mean
    correct = area_range.include?(true_ave) 
    
    print_data && print(correct ? " OK\n" : " NG\n")
    return correct
  end

  ## 累積分布からt分布信頼区間pの値を求める
  ## 挟み撃ち的にできるか？
  # fは自由度 nをサンプル数として区別する

  # 右側
  def t_inv_r(f, p, max: 10, step: 1e-4)
    x = max
    
    # 最初のwhileを通すよう初期値設定
    t = 1.1
    
    while t > p
      t = Distribution::T.cdf(x, f)
      x -= step
    end
    return x
  end

  # 左側
  def t_inv_l(f, p, min: -10, step: 1e-4)
    x = min
    
    # 最初のwhileを通すよう初期値設定
    t = -0.1
    
    while t < 1 - p
      t = Distribution::T.cdf(x, f)
      x += step
    end
    return x
  end

  # 両側
  def t_inv(f, p, max:10, step: 1e-4)
    x = max
    
    # 最初のwhileを通すよう初期値設定
    t = 1.1
    
    while t > 1.0 - (1.0 - p) / 2.0
      t = Distribution::T.cdf(x, f)
      x -= step
    end
    return [x, -x]
  end

  ## 累積分布からχ2分布信頼区間pの値を求める
  ## 挟み撃ち的にできるか？
  # fは自由度 nをサンプル数として区別する

  # 右側
  def chi_inv_r(f, p, max: 10, step: 1e-4)
    x = max
    
    # 最初のwhileを通すよう初期値設定
    chi = 1.1
    
    while chi > p
      chi = Distribution::ChiSquare.cdf(x, f)
      x -= step
    end
    return x
  end

  # 左側
  def chi_inv_l(f, p, min: 0, step: 1e-4)
    x = min
    
    # 最初のwhileを通すよう初期値設定
    chi = -0.1
    
    while chi < 1 - p
      chi = Distribution::ChiSquare.cdf(x, f)
      x += step
    end
    return x
  end

  # 両側
  def chi_inv(f, p, max:10, step: 1e-4)
    right = chi_inv_r(f, p, max: max, step: step)
    left = chi_inv_l(f, p, min: 0, step: step)
    return [right, left]
  end
end
