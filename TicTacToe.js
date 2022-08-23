// Tic Tace Toe Game
// Lekhraj Sharma
// Experiments in ReactJS
// Aug 2022

// square of game board
// pure component (doe snot derive from a React component)
function Square(props) {
  return (
    <button
      className="square"
      style={{ backgroundColor: props.bcolor }} // highlight color
      onClick={props.onClick}
    >
      {props.value}
    </button>
  );
}

// Game Board
// -------------
// | 0 | 1 | 2 |
// | 3 | 4 | 5 |
// | 6 | 7 | 8 |
// -------------

class Board extends React.Component {
  renderSquare(i) {
    // Helper to render an square
    // highlight background color
    const bcolor = (this.props.winnerCells.includes(i) ? "yellow" : "white");
    return (
      <Square
        value={this.props.squares[i]}
        onClick={() => this.props.onClick(i)}
        bcolor={bcolor}
      />
    );
  }

  render() {
    return (
      <div>
        <div className="board-row">
          {this.renderSquare(0)}
          {this.renderSquare(1)}
          {this.renderSquare(2)}
        </div>
        <div className="board-row">
          {this.renderSquare(3)}
          {this.renderSquare(4)}
          {this.renderSquare(5)}
        </div>
        <div className="board-row">
          {this.renderSquare(6)}
          {this.renderSquare(7)}
          {this.renderSquare(8)}
        </div>
      </div>
    );
  }
}

// Helper function to see if we have an winner!
// We have an winner in our 3x3 board, if any of rows or columns or diagonals have same mark!
// -------------
// | 0 | 1 | 2 |
// | 3 | 4 | 5 |
// | 6 | 7 | 8 |
// -------------
// If we do find an winner, we return both the winning player as well as winning squares

function calculateWinner(squares) {
  const lines = [
    [0, 1, 2], // top row
    [3, 4, 5], // middle row
    [6, 7, 8], // bottom row
    [0, 3, 6], // first column
    [1, 4, 7], // second column
    [2, 5, 8], // third column
    [0, 4, 8], // left diagonal
    [2, 4, 6] // right diagonal
  ];

  for (let i = 0; i < lines.length; i++) {
    const [a, b, c] = lines[i];
    if (squares[a] && squares[a] === squares[b] && squares[a] === squares[c]) {
      return [squares[a], lines[i]];
    }
  }
  return [null, [null, null, null]];
}

// Helper to remove [last+1..length] from an object array
// Normal splice method do not really work for state object arrays!
function removeElements(arr, last) {
  // Simulates arr.splice(last+1)
  for (let j = arr.length - 1; j > last; j--) {
    arr.pop();
  }
}

// Tic Tac Toe Game
// <Status Bar>
// <Board>
// <Move List to goback in history>

class Game extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      history: [
        {
          squares: Array(9).fill(null)
        }
      ],
      stepNumber: 0,
      xIsNext: true,
      winner: null,
      winnerCells: Array(3).fill(null)
    };
  }

  handleClick(i) {
    const history = this.state.history.slice();
    removeElements(history, this.state.stepNumber);
    const current = history[this.state.stepNumber];
    const squares = current.squares.slice();
    if (this.state.winner || squares[i]) {
      return;
    }
    squares[i] = this.state.xIsNext ? "X" : "O";
    this.setState({
      history: history.concat([
        {
          squares: squares
        }
      ]),
      stepNumber: history.length,
      xIsNext: !this.state.xIsNext
    });
  }

  // Helper to change state when we move back in history
  jumpTo(step) {
    this.setState({
      stepNumber: step,
      xIsNext: step % 2 == 0
    });
  }

  render() {
    const history = this.state.history.slice();
    removeElements(history, this.state.stepNumber);

    const current = history[this.state.stepNumber];

    [this.state.winner, this.state.winnerCells] = calculateWinner(
      current.squares
    );

    // Build drop down list item list
    const moves = history.map((step, move) => {
      const desc = move ? "Go to move #" + move : "Go to game start";
      return (
        // Drop down list item
        <option key={move} value={move}>
           {desc}
        </option>
      );
    });

    // Game status header
    let status;
    if (this.state.winner) {
      status = "Winner: " + this.state.winner;
    } else {
      status = "Next player: " + (this.state.xIsNext ? "X" : "O");
    }

    // <Status Bar>
    // <Board>
    // <Move Dropdown List to goback in history>
    return (
      <div className="game">
        <div className="game-board">
          <div>{status}</div>
          <Board
            squares={current.squares}
            onClick={(i) => this.handleClick(i)}
            winnerCells={this.state.winnerCells}
          />
        </div>
        <div className="game-info">
          <select onChange={(event) => this.jumpTo(event.target.value)}>
            {moves}
          </select>
        </div>
      </div>
    );
  }
}

// ========================================

const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<Game />); // start rendering our Game class
