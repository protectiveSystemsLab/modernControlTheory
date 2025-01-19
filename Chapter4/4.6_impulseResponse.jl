using LinearAlgebra
using Plots
Plots.gr()				# バックエンドの指定（デフォルトはGR）

# 単位系 t, m , kN
m = 1.0				    # 単位質量 1.0 ton
T = 1.0				    # 固有周期
w = 2*pi/T				# 固有円振動数
h = 0.05				# 減衰定数
k = (4*pi^2/T^2)*m			# ばね剛性
c = 2*h*w*m				# 減衰係数

A1 = [0 1; -k/m -c/m]	# 状態空間方程式の係数行列
b1 = [0; -1]
c1 = [1 0]

dt = 0.02				# 時間刻み
n = 501				    # 時刻データ数
t = dt*(0:n-1)			# 時刻配列

sol    = eigen(A1)		# 固有値解析
lambda = sol.values		# 固有値
U1     = sol.vectors	# 固有モード行列
U1inv  = inv(U1)

z1 = zeros(n,1)			# 変位データ配列

for j=1:n
    eAt = [exp(lambda[1]*t[j]) 0;
           0 exp(lambda[2]*t[j])]	# 状態遷移行列
    z1[j] = real(c1 * U1 * eAt * U1inv * b1)[1]
    # 上記は複素数型の解となるが虚部は消えるので，実部のみを
    # 取り出す．また，要素が一つの列ベクトルになるので
    # "[1]"によりスカラー値を取り出す．
end

# 時刻歴応答の図化
p=plot(t,z1,
    title="Time History",
    size=(600,300),
    label="displacement",
    xlabel="time (s)",
    ylabel="disp. (m)",
    xlims = (0,10),
    ylims = (-Inf, Inf),
    lw = 2)

display(p)